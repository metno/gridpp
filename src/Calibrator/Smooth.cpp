#include "Smooth.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
CalibratorSmooth::CalibratorSmooth(Variable::Type iVariable):
      Calibrator(),
      mSmoothRadius(3),
      mVariable(iVariable) {
}

bool CalibratorSmooth::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& precip = *iFile.getField(mVariable, t);
      Field precipRaw = precip;

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {

            for(int e = 0; e < nEns; e++) {
               float total = 0;
               int count = 0;
               for(int ii = std::max(0, i-mSmoothRadius); ii <= std::min(nLat-1, i+mSmoothRadius); ii++) {
                  for(int jj = std::max(0, j-mSmoothRadius); jj <= std::min(nLon-1, j+mSmoothRadius); jj++) {
                     float value = precipRaw(ii,jj,e);
                     if(Util::isValid(value)) {
                        total += value;
                        count++;
                     }
                  }
               }
               if(count > 0) {
                  precip(i,j,e) = total / count;
               }
               else {
                  precip(i,j,e) = Util::MV;
                  // std::stringstream ss;
                  // ss << "No precip value computed for [" << t << "," << i << "," << j << "," << e << "]";
                  // Util::warning(ss.str());
               }
            }
         }
      }
   }
   return true;
}

void CalibratorSmooth::setSmoothRadius(int iNumPoints) {
   if(!Util::isValid(iNumPoints) || iNumPoints < 1) {
      std::stringstream ss;
      ss << "CalibratorSmooth: Smoothing radius must be >= 1";
      Util::error(ss.str());
   }
   mSmoothRadius = iNumPoints;
}

int CalibratorSmooth::getSmoothRadius() const {
   return mSmoothRadius;
}
std::string CalibratorSmooth::description() {
   std::stringstream ss;
   ss << "   -c smooth                    Smooth the ensemble by averaging spatially, using a neighbourhood." << std::endl;
   ss << "      smoothRadius=3            Average gridpoints within a neighbourhood box of points within" << std::endl;
   ss << "                                +- radius in both east-west and north-south direction." << std::endl;
   return ss.str();
}
