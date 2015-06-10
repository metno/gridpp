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
      mOperator(Util::OperatorMean),
      mQuantile(Util::MV) {
   iOptions.getValue("radius", mRadius);
   if(mRadius < 0) {
      std::stringstream ss;
      ss << "CalibratorNeighbourhood: 'radius' (" << mRadius << ") must be >= 0";
      Util::error(ss.str());
   }

   std::string op;
   if(iOptions.getValue("operator", op)) {
      if(op == "mean") {
         mOperator = Util::OperatorMean;
      }
      else if(op == "min") {
         mOperator = Util::OperatorQuantile;
         mQuantile = 0;
      }
      else if(op == "max") {
         mOperator = Util::OperatorQuantile;
         mQuantile = 1;
      }
      else if(op == "median") {
         mOperator = Util::OperatorQuantile;
         mQuantile = 0.5;
      }
      else if(op == "std") {
         mOperator = Util::OperatorStd;
      }
      else if(op == "quantile"){
         mOperator = Util::OperatorQuantile;
         if(!iOptions.getValue("quantile", mQuantile)) {
            Util::error("CalibratorNeighbourhood: option 'quantile' is required");
         }
         if(!Util::isValid(mQuantile) || mQuantile < 0 || mQuantile > 1) {
            Util::error("CalibratorNeighbourhood: 'quantile' must be on the interval [0,1]");
         }
      }
      else {
         Util::error("CalibratorNeighbourhood: Unrecognized value for 'operator'");
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
               precip(i,j,e) = Util::applyOperator(neighbourhood, mOperator, mQuantile);
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
   ss << Util::formatDescription("-c neighbourhood", "Applies an operator on a neighbourhood (example by averaging across a neighbourhood thereby smoothing the field).") << std::endl;
   ss << Util::formatDescription("   radius=3", "Use gridpoints within this number of points within in both east-west and north-south direction.") << std::endl;
   ss << Util::formatDescription("   operator=mean", "What operator should be applied to the neighbourhood? One of 'mean', 'median', 'min', 'max', 'std', or 'quantile'. 'std' is the population standard deviation.") << std::endl;
   ss << Util::formatDescription("   quantile=undef", "If operator=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
   return ss.str();
}
