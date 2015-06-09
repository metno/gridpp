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
      mOperator(OperatorMean),
      mQuantile(Util::MV) {
   iOptions.getValue("radius", mRadius);
   if(mRadius < 0) {
      std::stringstream ss;
      ss << "CalibratorSmooth: Smoothing radius (" << mRadius << ") must be >= 1";
      Util::error(ss.str());
   }

   std::string op;
   if(iOptions.getValue("operator", op)) {
      if(op == "mean") {
         mOperator = OperatorMean;
      }
      else if(op == "min") {
         mOperator = OperatorQuantile;
         mQuantile = 0;
      }
      else if(op == "max") {
         mOperator = OperatorQuantile;
         mQuantile = 1;
      }
      else if(op == "median") {
         mOperator = OperatorQuantile;
         mQuantile = 0.5;
      }
      else if(op == "std") {
         mOperator = OperatorStd;
      }
      else if(op == "quantile"){
         mOperator = OperatorQuantile;
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
               precip(i,j,e) = compute(neighbourhood, mOperator, mQuantile);
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
float CalibratorNeighbourhood::compute(const std::vector<float>& neighbourhood, OperatorType iOperator, float iQuantile) {
   // Initialize to missing
   float value = Util::MV;
   if(iOperator == OperatorMean) {
      float total = 0;
      int count = 0;
      for(int n = 0; n < neighbourhood.size(); n++) {
         if(Util::isValid(neighbourhood[n])) {
            total += neighbourhood[n];
            count++;
         }
      }
      if(count > 0) {
         value = total / count;
      }
   }
   else if(iOperator == OperatorStd) {
      // STD = sqrt(E[X^2] - E[X]^2)
      // The above formula is unstable when the variance is small and the mean is large.
      // Use the property that VAR(X) = VAR(X-K). Provided K is any element in the neighbourhood
      // the resulting calculation of VAR(X-K) is stable. Set K to the first non-missing value.
      float total  = 0;
      float total2 = 0;
      float K = Util::MV;
      int count = 0;
      for(int n = 0; n < neighbourhood.size(); n++) {
         if(Util::isValid(neighbourhood[n])) {
            if(!Util::isValid(K))
               K = neighbourhood[n];
            assert(Util::isValid(K));

            total  += neighbourhood[n] - K;
            total2 += (neighbourhood[n] - K)*(neighbourhood[n] - K);
            count++;
         }
      }
      if(count > 0) {
         float mean  = total / count;
         float mean2 = total2 / count;
         float var   = mean2 - mean*mean;
         if(var < 0) {
            // This should never happen
            var = 0;
            Util::warning("CalibratorNeighbourhood: Problems computing std, unstable result. Setting value to 0");
         }
         float std = sqrt(var);
         value = std;
      }
   }
   else if(iOperator == OperatorQuantile) {
      // Remove missing
      std::vector<float> cleanHood;
      cleanHood.reserve(neighbourhood.size());
      for(int i = 0; i < neighbourhood.size(); i++) {
         if(Util::isValid(neighbourhood[i]))
            cleanHood.push_back(neighbourhood[i]);
      }
      int N = cleanHood.size();
      if(N > 0) {
         std::sort(cleanHood.begin(), cleanHood.end());
         int lowerIndex = floor(iQuantile * (N-1));
         int upperIndex = ceil(iQuantile * (N-1));
         float lowerQuantile = (float) lowerIndex / (N-1);
         float upperQuantile = (float) upperIndex / (N-1);
         float lowerValue = cleanHood[lowerIndex];
         float upperValue = cleanHood[upperIndex];
         if(lowerIndex == upperIndex) {
            value = lowerValue;
         }
         else {
            assert(upperQuantile > lowerQuantile);
            assert(iQuantile >= lowerQuantile);
            float f = (iQuantile - lowerQuantile)/(upperQuantile - lowerQuantile);
            assert(f >= 0);
            assert(f <= 1);
            value   = lowerValue + (upperValue - lowerValue) * f;
         }
      }
   }
   return value;
}
