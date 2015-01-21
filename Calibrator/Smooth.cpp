#include "Smooth.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
CalibratorSmooth::CalibratorSmooth(Variable::Type iMainPredictor):
      Calibrator(),
      mSmoothRadius(3),
      mMainPredictor(iMainPredictor) {
}

void CalibratorSmooth::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   int numInvalidRaw = 0;
   int numInvalidCal = 0;

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   DownscalerNearestNeighbour::getNearestNeighbourFast(iFile, iFile, nearestI, nearestJ);

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& precip = *iFile.getField(Variable::Precip, t);
      Field precipRaw = precip;

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int Icenter = nearestI[i][j];
            int Jcenter = nearestJ[i][j];

            for(int e = 0; e < nEns; e++) {
               float total = 0;
               int count = 0;
               for(int ii = std::max(0, Icenter-mSmoothRadius); ii <= std::min(nLat-1, Icenter+mSmoothRadius); ii++) {
                  for(int jj = std::max(0, Jcenter-mSmoothRadius); jj <= std::min(nLon-1, Jcenter+mSmoothRadius); jj++) {
                     float value = precipRaw[ii][jj][e];
                     if(Util::isValid(value)) {
                        total += value;
                        count++;
                     }
                  }
               }
               if(count > 0) {
                  precip[i][j][e] = total / count;
               }
               else {
                  precip[i][j][e] = Util::MV;
                  std::stringstream ss;
                  ss << "No precip value computed for [" << t << "," << i << "," << j << "," << e << "]";
                  Util::warning(ss.str());
               }
            }
         }
      }
   }
}

void CalibratorSmooth::setSmoothRadius(int iNumPoints) {
   mSmoothRadius = iNumPoints;
}

std::string CalibratorSmooth::description() {
   std::stringstream ss;
   ss << "   -c smooth                    Smooth the ensemble by averaging spatially, using a neighbourhood." << std::endl;
   ss << "      smoothRadius=3            Average gridpoints within a neighbourhood box of points within" << std::endl;
   ss << "                                +- radius in both east-west and north-south direction." << std::endl;
   return ss.str();
}
