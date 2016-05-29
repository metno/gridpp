#include "Coastal.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorCoastal::CalibratorCoastal(Variable::Type iVariable, const Options& iOptions) :
      Calibrator(iOptions),
      mVariable(iVariable),
      mSearchRadius(3) {
   iOptions.getValue("searchRadius", mSearchRadius);
}
bool CalibratorCoastal::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();
   vec2 lsm = iFile.getLandFractions();

   if(iParameterFile->getNumParameters() == 0) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must have at least one dataacolumns");
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Parameters parameters;
      if(!iParameterFile->isLocationDependent())
         parameters = iParameterFile->getParameters(t);
      const FieldPtr field = iFile.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            if(iParameterFile->isLocationDependent())
               parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
            for(int e = 0; e < nEns; e++) {
               float lowerValue = Util::MV;
               float upperValue = Util::MV;
               float minLsm = 1;
               float maxLsm = 0;
               float nearestNeighbour = (*field)(i,j,e);
               // Look inside a neighbourhood
               for(int ii = std::max(0, i-mSearchRadius); ii <= std::min(nLat-1, i+mSearchRadius); ii++) {
                  for(int jj = std::max(0, j-mSearchRadius); jj <= std::min(nLon-1, j+mSearchRadius); jj++) {
                     float currLsm = lsm[ii][jj];
                     if(currLsm < minLsm) {
                        lowerValue = (*field)(ii,jj,e);
                        minLsm = currLsm;
                     }
                     if(currLsm > minLsm) {
                        upperValue = (*field)(ii,jj,e);
                        maxLsm = currLsm;
                     }
                  }
               }

               if(Util::isValid(lowerValue) && Util::isValid(upperValue)) {
                  float a = parameters[0];
                  float b = parameters[0];
                  float c = parameters[0];
                  float range = upperValue - lowerValue;
                  float value = a + b * nearestNeighbour + c * range;
                  (*field)(i,j,e)  = value;
               }
               else {
                  (*field)(i,j,e) = nearestNeighbour;
               }
            }
         }
      }
   }
   return true;
}

Parameters CalibratorCoastal::train(const std::vector<ObsEnsField>& iData, int iIobs, int iJobs, int iIEns, int iJEns) const {
   /*
   float totalObs = 0;
   float totalForecast = 0;
   float totalForecast2 = 0;
   float totalObsForecast = 0;
   int counter = 0;
   // Compute predictors in model
   for(int i = 0; i < iData.size(); i++) {
      float obs = iData[i].first;
      std::vector<float> ens = iData[i].second;
      float ensMean = Util::calculateStat(ens, Util::StatTypeMean);
      if(Util::isValid(obs) && Util::isValid(ensMean)) {
         totalObs += obs;
         totalForecast += ensMean;
         totalForecast2 += ensMean*ensMean;
         totalObsForecast += obs*ensMean;
         counter++;
      }
   }

   if(counter <= 0) {
      std::stringstream ss;
      ss << "CalibratorCoastal: Cannot train regression, no valid data";
      Util::error(ss.str());
   }

   std::vector<float> values;
   float meanObs = totalObs / counter;
   float meanForecast = totalForecast / counter;
   float meanForecast2 = totalForecast2 / counter;
   float meanObsForecast = totalObsForecast / counter;
   if(mOrder == 0) {
      values.push_back(meanObs);
   }
   else if(mOrder == 1) {
      float intercept;
      float slope;
      if(mIntercept) {
         slope = (meanObsForecast - meanForecast*meanObs)/(meanForecast2 - meanForecast*meanForecast);
         intercept = meanObs - slope*meanForecast;
      }
      else {
         intercept = 0;
         slope = meanObsForecast / meanForecast2;
      }
      values.push_back(intercept);
      values.push_back(slope);
   }
   else {
      abort();
   }

   Parameters par(values);

   return par;
  */
}

std::string CalibratorCoastal::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c coastal", "Weights land and sea forecasts") << std::endl;
   ss << Util::formatDescription("   searchRadius=3", "") << std::endl;
   return ss.str();
}
