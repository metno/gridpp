#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"

// For each gridpoint in iTo, find the nearest neighbour from iFrom
// Output these into iI and iJ matrices, which represent the i and j indicies into the
// grid in iTo
void getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);
void getNearestNeighbourFast(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);

int main(int argc, const char *argv[]) {
   if(argc != 2 && argc != 3) {
      std::cout << "usage: nearestNeighbour file1 [file2]" << std::endl;
      return 1;
   }
   std::string inputFile     = argv[1];
   std::string outputFile    = argv[1];
   if(argc == 3)
      inputFile = argv[2];

   const FileArome ifile(inputFile);
   FileArome ofile(outputFile);

   std::cout << "Input size:  " << ifile.getNumLat() << "x" << ifile.getNumLon() << std::endl;
   std::cout << "Output size: " << ofile.getNumLat() << "x" << ofile.getNumLon() << std::endl;

   vec2Int I1, I2, J1, J2;
   double t0 = Util::clock();
   getNearestNeighbourFast(ifile, ofile, I2, J2);
   double t1 = Util::clock();
   std::cout << "Fast: " << t1-t0 << std::endl;
   double t2 = Util::clock();
   getNearestNeighbour(ifile, ofile, I1, J1);
   double t3 = Util::clock();
   std::cout << "Slow: " << t3-t2 << std::endl;
   float speedup = (t3-t2)/(t1-t0);
   std::cout << "Speedup factor: " << speedup << "x" << std::endl;


   assert(I1.size() == I2.size());
   assert(J1.size() == J2.size());
   for(int i = 0; i < I1.size(); i++) {
      for(int j = 0; j < I1[0].size(); j++) {
         if(I1[i][j] != I2[i][j] || J1[i][j] != J2[i][j]) {
            std::cout << "ERROR: Different results for (" << i << " " << j << "): " << I1[i][j] << " " << J1[i][j] << " " << I2[i][j] << " " << J2[i][j] << std::endl;
         }
      }
   }
}

void getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumLon();
   int nLat = iTo.getNumLat();

   iI.resize(nLat);
   iJ.resize(nLat);

   // Loop over points in the output grid (i and j)
   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      iI[i].resize(nLon, Util::MV);
      iJ[i].resize(nLon, Util::MV);
      for(int j = 0; j < nLon; j++) {
         if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
            // Find the nearest neighbour from input grid (ii, jj)
            float minDist = Util::MV;
            for(int ii = 0; ii < ilats.size(); ii++) {
               for(int jj = 0; jj < ilats[0].size(); jj++) {
                  if(Util::isValid(ilats[ii][jj]) && Util::isValid(ilons[ii][jj])) {
                     float currDist = Util::getDistance(ilats[ii][jj], ilons[ii][jj], olats[i][j], olons[i][j]);
                     if(!Util::isValid(minDist) || currDist < minDist) {
                        iI[i][j]= ii;
                        iJ[i][j]= jj;
                        minDist = currDist;
                     }
                  }
               }
            }
         }
      }
   }
}

void getNearestNeighbourFast(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumLon();
   int nLat = iTo.getNumLat();

   iI.resize(nLat);
   iJ.resize(nLat);

   float tol = 0.2;

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      int I = 0;
      int J = 0;
      iI[i].resize(nLon, 0);
      iJ[i].resize(nLon, 0);
      for(int j = 0; j < nLon; j++) {
         int counter = 0;
         float currLat = olats[i][j];
         float currLon = olons[i][j];
         while(true) {
            if(fabs(ilons[I][J] - currLon) < tol && fabs(ilats[I][J] - currLat) < tol) {
               int startI = std::max(0, I-10);
               int startJ = std::max(0, J-10);
               int endI   = std::min(iFrom.getNumLat()-1, I+10);
               int endJ   = std::min(iFrom.getNumLon()-1, J+10);
               float minDist = Util::MV;
               for(int ii = startI; ii <= endI; ii++) {
                  for(int jj = startJ; jj <= endJ; jj++) {
                     float currDist = Util::getDistance(olats[i][j], olons[i][j], ilats[ii][jj], ilons[ii][jj]);
                     if(!Util::isValid(minDist) || currDist < minDist) {
                        I = ii;
                        J = jj;
                        minDist = currDist;
                     }
                  }
               }
               break;
            }
            else {
               assert(I >= 0);
               assert(J >= 0);
               assert(ilons.size() > I);
               assert(ilons[I].size() > J);
               assert(olons.size() > i);
               assert(olons[i].size() > j);
               assert(ilats.size() > I);
               assert(ilats[I].size() > J);
               assert(olats.size() > i);
               assert(olats[i].size() > j);
               if(ilons[I][J] < currLon-tol)
                  J++;
               else if(ilons[I][J] > currLon+tol)
                  J--;
               if(ilats[I][J] < currLat-tol)
                  I++;
               else if(ilats[I][J] < currLat+tol)
                  I--;
               I = std::min(iFrom.getNumLat()-1, std::max(0, I));
               J = std::min(iFrom.getNumLon()-1, std::max(0, J));
            }
            counter++;
            if(counter > 1000) {
               float minDist = Util::MV;
               for(int ii = 0; ii < iFrom.getNumLat(); ii++) {
                  for(int jj = 0; jj < iFrom.getNumLon(); jj++) {
                     float currDist = Util::getDistance(olats[i][j], olons[i][j], ilats[ii][jj], ilons[ii][jj]);
                     if(Util::isValid(currDist) && (!Util::isValid(minDist) || currDist < minDist)) {
                        I = ii;
                        J = jj;
                        minDist = currDist;
                     }
                  }
               }
               break;
            }
         }
         iI[i][j] = I;
         iJ[i][j] = J;
      }
   }
}
