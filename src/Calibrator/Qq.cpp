#include "Qq.h"
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
#include "gridpp.h"
#include <cmath>

CalibratorQq::CalibratorQq(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mLowerQuantile(0),
      mUpperQuantile(1),
      mX(Util::MV),
      mY(Util::MV),
      mPolicy(gridpp::OneToOne) {

   std::string extrapolation;
   if(iOptions.getValue("extrapolation", extrapolation)) {
      if(extrapolation == "1to1")
         mPolicy = gridpp::OneToOne;
      else if(extrapolation == "meanSlope")
         mPolicy = gridpp::MeanSlope;
      else if(extrapolation == "nearestSlope")
         mPolicy = gridpp::NearestSlope;
      else if(extrapolation == "zero")
         mPolicy = gridpp::Zero;
      else {
         Util::error("CalibratorQq: value for 'extrapolation' not recognized");
      }
   }
   iOptions.getValues("quantiles", mQuantiles);
   iOptions.getValues("extraObs", mExtraObs);
   iOptions.getValues("extraFcst", mExtraFcst);
   if(mExtraObs.size() != mExtraFcst.size()) {
      Util::error("extraObs must have the same length as extraFcst");
   }
   iOptions.getValue("x", mX);
   iOptions.getValue("y", mY);
   iOptions.check();
}
bool CalibratorQq::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(iParameterFile->getNumParameters() % 2 != 0) {
      std::stringstream ss;
      ss << "Parameter file '" + iParameterFile->getFilename() + "' must have an even number of datacolumns (currently " << iParameterFile->getNumParameters() << " parameters)";
      Util::error(ss.str());
   }

   const int nEns = iFile.getNumEns();
   const int nTime = iFile.getNumTime();
   gridpp::Grid grid(iFile.getLats(), iFile.getLons());

   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      gridpp::Parameters parameters = iParameterFile->getApiParameters(t);
      for(int e = 0; e < nEns; e++) {
         field->set(gridpp::quantile_mapping(grid, (*field)(e), mPolicy, parameters), e);
      }
   }
   return true;
}

Parameters CalibratorQq::train(const std::vector<ObsEns>& iData) const {
   double timeStart = Util::clock();
   std::vector<float> obs, fcst;
   std::vector<float> values;
   values.reserve(2*iData.size());
   int counter = 0;
   // Compute predictors in model
   for(int i = 0; i < iData.size(); i++) {
      float currObs = iData[i].first;
      std::vector<float> ens = iData[i].second;
      float currFcst = Util::calculateStat(ens, Util::StatTypeMean);
      if(Util::isValid(currObs) && Util::isValid(currFcst)) {
         obs.push_back(currObs);
         fcst.push_back(currFcst);
      }
   }
   assert(obs.size() == fcst.size());

   if(obs.size() == 0 || fcst.size() == 0) {
      // std::stringstream ss;
      // ss << "CalibratorQq: No valid data, no correction will be made.";
      // Util::warning(ss.str());
   }

   // Sort
   std::sort(obs.begin(), obs.end());
   std::sort(fcst.begin(), fcst.end());
   std::vector<float> obsQ;
   std::vector<float> fcstQ;
   if(obs.size() == 0) {

   }
   else if(mQuantiles.size() > 0) {
      int N = obs.size();
      std::vector<float> x(N);
      for(int i = 0; i < N; i++) {
         x[i] = (float) i / (N-1);
      }
      for(int i = 0; i < mQuantiles.size(); i++) {
         float currObs  = Util::interpolate(mQuantiles[i], x, obs);
         float currFcst = Util::interpolate(mQuantiles[i], x, fcst);
         obsQ.push_back(currObs);
         fcstQ.push_back(currFcst);
      }
   }
   else {
      for(int i = 0; i < obs.size(); i++) {
         obsQ.push_back(obs[i]);
         fcstQ.push_back(fcst[i]);
      }
   }

   // Add extra points
   if(mExtraObs.size() > 0) {
      for(int i = 0; i < mExtraObs.size(); i++) {
         obsQ.push_back(mExtraObs[i]);
         fcstQ.push_back(mExtraFcst[i]);
      }
      std::sort(obsQ.begin(), obsQ.end());
      std::sort(fcstQ.begin(), fcstQ.end());
   }

   // Put into parameters
   for(int i = 0; i < obsQ.size(); i++) {
      values.push_back(obsQ[i]);
      values.push_back(fcstQ[i]);
   }

   Parameters par(values);

   double timeEnd = Util::clock();
   // std::cout << "Time: " << timeEnd - timeStart << std::endl;
   return par;
}

std::string CalibratorQq::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c qq", "Quantile-quantile mapping. Calibrates forecasts based on a map of sorted observations and sorted forecasts. For a given raw forecast, the quantile within the historical forecasts is found. Then the observation at the same quantile is used as the calibrated forecast. A parameter file is required with an even number of columns as follows") << std::endl;
      ss << Util::formatDescription("", "[obs0 fcs0 obs1 fcst1 .. obsN fcstN") << std::endl;
      ss << Util::formatDescription("", "Note that observations and forecasts must be sorted. I.e obs0 does not necessarily correspond to the time when fcst0 was issued. If a parameter set has one or more missing values, then the ensemble using this parameter set is not processed.") << std::endl;
      ss << Util::formatDescription("   extrapolation=1to1", "If a forecast is outside the curve, how should extrapolation be done? '1to1': Use a slope of 1, i.e. preserving the bias at the nearest end point; 'meanSlope': Use the average slope from the lower to upper endpoints; 'nearestSlope': Use the slope through the two nearest points; 'zero': Use a slope of 0, meaning that the forecast will equal the max/min observation.") << std::endl;
      ss << Util::formatDescription("   quantiles=undef", "If creating training parameters, should specific quantiles be stored, instead of all observations and forecasts? Can be a vector of values between 0 and 1.") << std::endl;
      ss << Util::formatDescription("   extraObs=undef", "Add these extra observations points to the curve. Must be the same length as extraFcst.") << std::endl;
      ss << Util::formatDescription("   extraFcst=undef", "Add these extra forecasts points to the curve at the end. Must be the same length as extraObs.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c qq", "Quantile-quantile mapping") << std::endl;
   return ss.str();
}

void CalibratorQq::separate(const Parameters& iParameters, std::vector<float>& iObs, std::vector<float>& iFcst) const {
   iObs.clear();
   iFcst.clear();
   int N = iParameters.size() / 2;
   iObs.resize(N, Util::MV);
   iFcst.resize(N, Util::MV);
   for(int i = 0; i < N; i++) {
      iObs[i] = iParameters[2*i];
      iFcst[i] = iParameters[2*i+1];
   }
   if(mExtraObs.size() > 0) {
      for(int i = 0; i < mExtraObs.size(); i++) {
         iObs.push_back(mExtraObs[i]);
         iFcst.push_back(mExtraFcst[i]);
      }
      std::sort(iObs.begin(), iObs.end());
      std::sort(iFcst.begin(), iFcst.end());
   }
}
