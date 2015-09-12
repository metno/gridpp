#include "Qq.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorQq::CalibratorQq(Variable::Type iVariable, const Options& iOptions) :
      Calibrator(iOptions),
      mVariable(iVariable),
      mLowerQuantile(0),
      mUpperQuantile(1),
      mPolicy(ExtrapolationPolicy::OneToOne) {

   std::string extrapolation;
   if(iOptions.getValue("extrapolation", extrapolation)) {
      if(extrapolation == "1to1")
         mPolicy = ExtrapolationPolicy::OneToOne;
      else if(extrapolation == "meanSlope")
         mPolicy = ExtrapolationPolicy::MeanSlope;
      else if(extrapolation == "nearestSlope")
         mPolicy = ExtrapolationPolicy::NearestSlope;
      else if(extrapolation == "zero")
         mPolicy = ExtrapolationPolicy::Zero;
      else {
         Util::error("CalibratorQq: value for 'extrapolation' not recognized");
      }
   }
}
bool CalibratorQq::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   const int nLat = iFile.getNumLat();
   const int nLon = iFile.getNumLon();
   const int nEns = iFile.getNumEns();
   const int nTime = iFile.getNumTime();
   const vec2 lats = iFile.getLats();
   const vec2 lons = iFile.getLons();
   const vec2 elevs = iFile.getElevs();

   if(iParameterFile->getNumParameters() % 2 != 0) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must have an even number of datacolumns");
   }

   for(int t = 0; t < nTime; t++) {
      // Retrieve the calibration parameters for this time
      // Overwrite them later if they are location dependent
      std::vector<float> obsVec, fcstVec;
      Parameters parameters;
      if(!iParameterFile->isLocationDependent()) {
         parameters = iParameterFile->getParameters(t);
         separate(parameters, obsVec, fcstVec);
      }
      const FieldPtr field = iFile.getField(mVariable, t);

      #pragma omp parallel for private(obsVec,fcstVec,parameters)
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int N = obsVec.size();
            if(iParameterFile->isLocationDependent()) {
               parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
               separate(parameters, obsVec, fcstVec);
            }
            // Only process if all parameters are valid
            bool hasValidParameters = true;
            for(int p = 0; p < parameters.size(); p++) {
               hasValidParameters = hasValidParameters && Util::isValid(parameters[p]);
            }
            if(hasValidParameters) {
               for(int e = 0; e < nEns; e++) {
                  float raw = (*field)(i,j,e);
                  float value = Util::MV;
                  if(Util::isValid(raw)) {
                     float smallestObs  = obsVec[0];
                     float smallestFcst = fcstVec[0];
                     float largestObs   = obsVec[obsVec.size()-1];
                     float largestFcst  = fcstVec[fcstVec.size()-1];

                     // Linear interpolation within curve
                     if(raw > smallestFcst && raw < largestFcst) {
                        value = Util::interpolate(raw, fcstVec, obsVec);
                     }
                     // Extrapolate outside curve
                     else {
                        float nearestObs;
                        float nearestFcst;
                        if(raw <= smallestFcst) {
                           nearestObs  = smallestObs;
                           nearestFcst = smallestFcst;
                        }
                        else {
                           nearestObs  = largestObs;
                           nearestFcst = largestFcst;
                        }
                        float slope = 1;
                        if(mPolicy == ExtrapolationPolicy::Zero) {
                           slope = 0;
                        }
                        if(mPolicy == ExtrapolationPolicy::OneToOne || N <= 1) {
                           slope = 1;
                        }
                        else if(mPolicy == ExtrapolationPolicy::MeanSlope) {
                           float dObs  = largestObs - smallestObs;
                           float dFcst = largestFcst - smallestFcst;
                           slope = dObs / dFcst;
                        }
                        else if(mPolicy == ExtrapolationPolicy::NearestSlope) {
                           float dObs;
                           float dFcst;
                           if(raw <= smallestFcst) {
                              dObs  = obsVec[1] - obsVec[0];
                              dFcst = fcstVec[1] - fcstVec[0];
                           }
                           else {
                              dObs  = obsVec[N-1] - obsVec[N-2];
                              dFcst = fcstVec[N-1] - fcstVec[N-2];
                           }
                           slope = dObs / dFcst;
                        }
                        value = nearestObs + slope * (raw - nearestFcst);
                     }
                     (*field)(i,j,e) = value;
                  }
                  else {
                     (*field)(i,j,e)  = Util::MV;
                  }
               }
            }
         }
      }
   }
   return true;
}

Parameters CalibratorQq::train(const TrainingData& iData, int iOffset) const {
   double timeStart = Util::clock();
   std::vector<ObsEns> data = iData.getData(iOffset);
   std::vector<float> obs, mean;
   std::vector<float> values;
   values.reserve(2*data.size());
   int counter = 0;
   // Compute predictors in model
   for(int i = 0; i < data.size(); i++) {
      float obs = data[i].first;
      std::vector<float> ens = data[i].second;
      float mean = Util::calculateStat(ens, Util::StatTypeMean);
      if(Util::isValid(obs) && Util::isValid(mean)) {
         values.push_back(obs);
         values.push_back(mean);
         counter++;
      }
   }

   if(counter <= 0) {
      std::stringstream ss;
      ss << "CalibratorQq: No valid data, no correction will be made.";
      Util::warning(ss.str());
   }

   Parameters par(values);

   double timeEnd = Util::clock();
   std::cout << "Time: " << timeEnd - timeStart << std::endl;
   return par;
}

std::string CalibratorQq::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c qq", "Quantile-quantile mapping. Calibrates forecasts based on a map of sorted observations and sorted forecasts. For a given raw forecast, the quantile within the historical forecasts is found. Then the observation at the same quantile is used as the calibrated forecast. A parameter file is required with an even number of columns as follows") << std::endl;
   ss << Util::formatDescription("", "[obs0 fcs0 obs1 fcst1 .. obsN fcstN") << std::endl;
   ss << Util::formatDescription("", "Note that observations and forecasts must be sorted. I.e obs0 does not necessarily correspond to the time when fcst0 was issued. If a parameter set has one or more missing values, then the ensemble using this parameter set is not processed.") << std::endl;
   ss << Util::formatDescription("   extrapolation=1to1", "If a forecast is outside the curve, how should extrapolation be done? '1to1': Use a slope of 1, i.e. preserving the bias at the nearest end point; 'meanSlope': Use the average slope from the lower to upper endpoints; 'nearestSlope': Use the slope through the two nearest points; 'zero': Use a slope of 0, meaning that the forecast will equal the max/min observation.") << std::endl;
   return ss.str();
}

void CalibratorQq::separate(const Parameters& iParameters, std::vector<float>& iObs, std::vector<float>& iFcst) {
   iObs.clear();
   iFcst.clear();
   int N = iParameters.size() / 2;
   iObs.resize(N, Util::MV);
   iFcst.resize(N, Util::MV);
   for(int i = 0; i < N; i++) {
      iObs[i] = iParameters[2*i];
      iFcst[i] = iParameters[2*i+1];
   }
}
