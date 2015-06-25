#include "Neighbourhood.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
CalibratorNeighbourhood::CalibratorNeighbourhood(Variable::Type iVariable, const Options& iOptions):
      Calibrator(),
      mRadius(3),
      mVariable(iVariable),
      mStatType(Util::StatTypeMean),
      mQuantile(Util::MV) {
   iOptions.getValue("radius", mRadius);
   if(mRadius < 0) {
      std::stringstream ss;
      ss << "CalibratorNeighbourhood: 'radius' (" << mRadius << ") must be >= 0";
      Util::error(ss.str());
   }

   std::string op;
   if(iOptions.getValue("stat", op)) {
      if(op == "mean") {
         mStatType = Util::StatTypeMean;
      }
      else if(op == "min") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 0;
      }
      else if(op == "max") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 1;
      }
      else if(op == "median") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 0.5;
      }
      else if(op == "std") {
         mStatType = Util::StatTypeStd;
      }
      else if(op == "quantile"){
         mStatType = Util::StatTypeQuantile;
         if(!iOptions.getValue("quantile", mQuantile)) {
            Util::error("CalibratorNeighbourhood: option 'quantile' is required");
         }
         if(!Util::isValid(mQuantile) || mQuantile < 0 || mQuantile > 1) {
            Util::error("CalibratorNeighbourhood: 'quantile' must be on the interval [0,1]");
         }
      }
      else {
         Util::error("CalibratorNeighbourhood: Unrecognized value for 'stat'");
      }
   }
}

bool CalibratorNeighbourhood::calibrateCore(File& iFile) const {
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
               // Put neighbourhood into vector
               std::vector<float> neighbourhood;
               int Ni = std::min(nLat-1, i+mRadius) - std::max(0, i-mRadius) + 1;
               int Nj = std::min(nLon-1, j+mRadius) - std::max(0, j-mRadius) + 1;
               assert(Ni > 0);
               assert(Nj > 0);
               neighbourhood.resize(Ni*Nj, Util::MV);
               int index = 0;
               for(int ii = std::max(0, i-mRadius); ii <= std::min(nLat-1, i+mRadius); ii++) {
                  for(int jj = std::max(0, j-mRadius); jj <= std::min(nLon-1, j+mRadius); jj++) {
                     float value = precipRaw(ii,jj,e);
                     assert(index < Ni*Nj);
                     neighbourhood[index] = value;
                     index++;
                  }
               }
               assert(index == Ni*Nj);
               precip(i,j,e) = Util::calculateStat(neighbourhood, mStatType, mQuantile);
            }
         }
      }
   }
   return true;
}

int CalibratorNeighbourhood::getRadius() const {
   return mRadius;
}

std::string CalibratorNeighbourhood::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c neighbourhood", "Applies a statistical operator on a neighbourhood (example by averaging across a neighbourhood thereby smoothing the field).") << std::endl;
   ss << Util::formatDescription("   radius=3", "Use gridpoints within this number of points within in both east-west and north-south direction.") << std::endl;
   ss << Util::formatDescription("   stat=mean", "What statistical operator should be applied to the neighbourhood? One of 'mean', 'median', 'min', 'max', 'std', or 'quantile'. 'std' is the population standard deviation.") << std::endl;
   ss << Util::formatDescription("   quantile=undef", "If stat=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
   return ss.str();
}
