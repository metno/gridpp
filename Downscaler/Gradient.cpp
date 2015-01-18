#include "Gradient.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerGradient::DownscalerGradient(Variable::Type iVariable) :
      Downscaler(),
      mRadius(3),
      mConstGradient(Util::MV),
      mVariable(iVariable) {
}

void DownscalerGradient::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumLat();
   int nLon = iOutput.getNumLon();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   vec2 ilats  = iInput.getLats();
   vec2 ilons  = iInput.getLons();
   vec2 ielevs = iInput.getElevs();
   vec2 olats  = iOutput.getLats();
   vec2 olons  = iOutput.getLons();
   vec2 oelevs = iOutput.getElevs();

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   DownscalerNearestNeighbour::getNearestNeighbourFast(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mVariable, t);
      Field& ofield = *iOutput.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int Icenter = nearestI[i][j];
            int Jcenter = nearestJ[i][j];
            for(int e = 0; e < nEns; e++) {
               float gradient = 0;
               if(!Util::isValid(mConstGradient)) {
                  // Compute gradient from neighbourhood
                  float meanXY  = 0; // elev*T
                  float meanX   = 0; // elev
                  float meanY   = 0; // T
                  float meanXX  = 0; // elev*elev
                  int   counter = 0;
                  for(int ii = std::max(0, Icenter-mRadius); ii <= std::min(iInput.getNumLat()-1, Icenter+mRadius); ii++) {
                     for(int jj = std::max(0, Jcenter-mRadius); jj <= std::min(iInput.getNumLon()-1, Jcenter+mRadius); jj++) {
                        float x = ielevs[ii][jj];
                        float y = ifield[ii][jj][e];
                        if(Util::isValid(x) && Util::isValid(y)) {
                           meanXY += x*y;
                           meanX  += x;
                           meanY  += y;
                           meanXX += x*x;
                           counter++;
                        }
                     }
                  }

                  if(fabs(meanXX - meanX*meanX) != 1) {
                     // Estimate lapse rate
                     meanXY /= counter;
                     meanX  /= counter;
                     meanY  /= counter;
                     meanXX /= counter;
                     gradient = -(meanXY - meanX*meanY)/(meanXX - meanX*meanX);
                  }
               }
               else {
                  gradient = mConstGradient;
               }

               float currElev = oelevs[i][j];
               float nearestElev = ielevs[Icenter][Jcenter];
               float dElev = currElev - nearestElev;
               // if(t == 0)
               //    std::cout << i << " " << j << " " << gradient << std::endl;

               ofield[i][j][e] = ifield[Icenter][Jcenter][e] + dElev * gradient;
            }
         }
      }
   }
}
void DownscalerGradient::setConstantGradient(float iGradient) {
   mConstGradient = iGradient;
}
void DownscalerGradient::setNeighbourhoodRadius(int iNumPoints) {
   mRadius = iNumPoints;
}
