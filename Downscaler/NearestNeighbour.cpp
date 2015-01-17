#include "NearestNeighbour.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerNearestNeighbour::DownscalerNearestNeighbour(Variable::Type iVariable) :
      Downscaler(),
      mVariable(iVariable) {
}

void DownscalerNearestNeighbour::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumLat();
   int nLon = iOutput.getNumLon();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   getNearestNeighbourFast(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mVariable, t);
      Field& ofield = *iOutput.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int I = nearestI[i][j];
            int J = nearestJ[i][j];
            for(int e = 0; e < nEns; e++) {
               ofield[i][j][e] = ifield[I][J][e];
            }
         }
      }
   }
}

void DownscalerNearestNeighbour::getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumLon();
   int nLat = iTo.getNumLat();

   iI.resize(nLat);
   iJ.resize(nLat);
   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      iI[i].resize(nLon, 0);
      iJ[i].resize(nLon, 0);
      for(int j = 0; j < nLon; j++) {
         if(j % 10 == 0)
            std::cout << i << " " << j << std::endl;
         float minDist = Util::MV;
         int I = 0;
         int J = 0;
         for(int ii = 0; ii < iFrom.getNumLat(); ii++) {
            for(int jj = 0; jj < iFrom.getNumLon(); jj++) {
               float currDist = Util::getDistance(olats[i][j], olons[i][j], ilats[ii][jj], ilons[ii][jj]);
               if(!Util::isValid(minDist) || currDist < minDist) {
                  I = ii;
                  J = jj;
                  minDist = currDist;
               }
            }
         }
         iI[i][j] = I;
         iJ[i][j] = J;
         assert(I < 100000);
         assert(J < 100000);
      }
   }
}

void DownscalerNearestNeighbour::getNearestNeighbourFast(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumLon();
   int nLat = iTo.getNumLat();

   iI.resize(nLat);
   iJ.resize(nLat);
   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      int I = 0;
      int J = 0;
      iI[i].resize(nLon, 0);
      iJ[i].resize(nLon, 0);
      for(int j = 0; j < nLon; j++) {
         if(j % 10 == 0)
            std::cout << i << " " << j << std::endl;
         int counter = 0;
         while(true) {
            // std::cout << "   " << ilats[I][J] << " " << ilons[I][J] << std::endl;
            if(fabs(ilons[I][J] - olons[i][j]) < 0.02 && fabs(ilats[I][J] - olats[i][j]) < 0.02) {
               int startI = std::max(0, I-10);
               int startJ = std::max(0, J-10);
               int endI   = std::min(iFrom.getNumLat()-1, I+10);
               int endJ   = std::min(iFrom.getNumLon()-1, J+10);
               // std::cout << i << " " << j << " " << olats[i][j] << " " << ilats[I][J] << " " << olons[i][j] << " " << ilons[I][J] << std::endl;
               // std::cout << i << " " << j << " " << abs(ilons[I][J] - olons[i][j]) << " " << abs(ilats[I][J] - olats[i][j]) << std::endl;
               // std::cout << i << " " << j << " Searching in [" << startI << " " << startJ << " " << endI << " " << endJ << "]" << std::endl;
               float minDist = Util::MV;
               for(int ii = startI; ii <= endI; ii++) {
                  for(int jj = startJ; jj <= endJ; jj++) {
                     float currDist = Util::getDistance(olats[i][j], olons[i][j], ilats[ii][jj], ilons[ii][jj]);
                     if(!Util::isValid(minDist) || currDist < minDist) {
                        // std::cout << ilats[ii][jj] << " " << ilons[ii][jj] << "    " << currDist << std::endl;
                        I = ii;
                        J = jj;
                        minDist = currDist;
                     }
                  }
               }
               // std::cout << "Found: " << i << " " << j << " " << olats[i][j] << " " << olons[i][j] << " " << ilats[I][J] << " " << ilons[I][J] << " " << minDist << std::endl;
               break;
            }
            else {
               if(ilons[I][J] < olons[i][j])
                  J++;
               else
                  J--;
               if(ilats[I][J] < olats[i][j])
                  I++;
               else
                  I--;
               I = std::min(iFrom.getNumLat()-1, std::max(0, I));
               J = std::min(iFrom.getNumLon()-1, std::max(0, J));
            }
            counter++;
            if(counter > 1000) {
               // std::cout << "Couldn't find for " << i << " " << j << std::endl;
               float minDist = Util::MV;
               for(int ii = 0; ii < iFrom.getNumLat(); ii++) {
                  for(int jj = 0; jj < iFrom.getNumLon(); jj++) {
                     float currDist = Util::getDistance(olats[i][j], olons[i][j], ilats[ii][jj], ilons[ii][jj]);
                     // std::cout << ii << " " << jj << " " << currDist << " " << minDist << std::endl;
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
         // std::cout << "Found: " << i << " " << j << " " << olats[i][j] << " " << olons[i][j] << " " << ilats[I][J] << " " << ilons[i][J] << std::endl;
      }
   }
}
