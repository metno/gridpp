#include "Neighbourhood.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "gridpp.h"
CalibratorNeighbourhood::CalibratorNeighbourhood(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mRadius(3),
      mStatType(Util::StatTypeMean),
      mStatTypeName("mean"),
      mFast(true),
      mNumThresholds(0),
      mQuantile(Util::MV) {
   iOptions.getValue("radius", mRadius);
   iOptions.getValue("fast", mFast);
   if(mRadius < 0) {
      std::stringstream ss;
      ss << "CalibratorNeighbourhood: 'radius' (" << mRadius << ") must be >= 0";
      Util::error(ss.str());
   }

   std::string op;
   if(iOptions.getValue("stat", mStatTypeName)) {
      bool status = Util::getStatType(mStatTypeName, mStatType);
      if(!status) {
         std::stringstream ss;
         ss << "Could not recognize stat=" << mStatTypeName;
         Util::error(ss.str());
      }
   }
   if(mStatType == Util::StatTypeQuantile) {
      iOptions.getValue("numThresholds", mNumThresholds);
      iOptions.getRequiredValue("quantile", mQuantile);
      if(!Util::isValid(mQuantile) || mQuantile < 0 || mQuantile > 1) {
         Util::error("'quantile' must be on the interval [0,1]");
      }
   }
   iOptions.check();
}

bool CalibratorNeighbourhood::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nTime = iFile.getNumTime();
   for(int t = 0; t < nTime; t++) {
      Field& output = *iFile.getField(mVariable, t);
      Field raw = output;
      if(iParameterFile != NULL) {
         if(iParameterFile->isLocationDependent()) {
            Util::error("Cannot use a location dependent parameter file for CalibratorNeighbourhood");
         }
         Parameters parameters = iParameterFile->getParameters(t);
         calibrateField(raw, output, &parameters);
      }
      else {
         calibrateField(raw, output);
      }
   }
   return true;
}
void CalibratorNeighbourhood::calibrateField(const Field& iInput, Field& iOutput, const Parameters* iParameters) const {
   double start_time = Util::clock();
   int radius = mRadius;
   int nEns = iInput.getNumEns();
   if(iParameters != NULL) {
      radius = (*iParameters)[0];
   }

   int count_stat = 0;
   for(int e = 0; e < nEns; e++) {
      const vec2 input = iInput(e);
      if(mStatType != Util::StatTypeQuantile && mStatType != Util::StatTypeMedian) {
          iOutput.set(gridpp::neighbourhood(input, radius, mStatTypeName), e);
      }
      else {
         iOutput.set(gridpp::neighbourhood_quantile(input, radius, mQuantile, mNumThresholds), e);
      }
   }
}

std::string CalibratorNeighbourhood::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c neighbourhood", "Applies a statistical operator on a neighbourhood (example by averaging across a neighbourhood thereby smoothing the field).") << std::endl;
      ss << Util::formatDescription("   radius=3", "Use gridpoints within this number of points within in both east-west and north-south direction. The radius can alternatively be specified using a location-independent parameter file, with one parameter.") << std::endl;
      ss << Util::formatDescription("   stat=mean", "What statistical operator should be applied to the neighbourhood? One of 'mean', 'median', 'min', 'max', 'quantile', 'std', or 'sum'. 'std' is the population standard deviation.") << std::endl;
      ss << Util::formatDescription("   quantile=undef", "If stat=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
      ss << Util::formatDescription("   numThresholds=0", "If stat=quantile, how many approximation thresholds to use? Use 0 for no approximation.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c neighbourhood", "Applies a statistical operator on a neighbourhood") << std::endl;
   return ss.str();
}
