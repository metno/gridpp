#include "Regression.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
#include "../TrainingData.h"
CalibratorRegression::CalibratorRegression(Variable::Type iVariable, const Options& iOptions) :
      Calibrator(iOptions),
      mVariable(iVariable),
      mOrder(1),
      mIntercept(true) {
   iOptions.getValue("order", mOrder);
   iOptions.getValue("intercept", mIntercept);
}
bool CalibratorRegression::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(iParameterFile == NULL) {
      Util::error("Calibrator 'regression' requires a valid parameter file");
   }
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   if(iParameterFile->getNumParameters() == 0) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must have at least one datacolumns");
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
               if(Util::isValid((*field)(i,j,e))) {
                  float total = 0;
                  // Accumulate a + b * fcst + c * fcst^2 ...
                  for(int p = 0; p < parameters.size(); p++) {
                     float coeff = parameters[p];
                     if(!Util::isValid(coeff)) {
                        total = Util::MV;
                        break;
                     }
                     total += coeff*pow((*field)(i,j,e), p);
                  }
                  (*field)(i,j,e)  = total;
               }
               else {
                  (*field)(i,j,e)  = Util::MV;
               }
            }
         }
      }
   }
   return true;
}

Parameters CalibratorRegression::train(const TrainingData& iData, int iOffset) const {
   if(mOrder > 1) {
      std::stringstream ss;
      ss << "CalibratorRegression: Cannot train regression of order greater than 1 (i.e a + bx)";
      Util::error(ss.str());
   }

   std::vector<ObsEns> data = iData.getData(iOffset);
   float totalObs = 0;
   float totalForecast = 0;
   float totalForecast2 = 0;
   float totalObsForecast = 0;
   int counter = 0;
   // Compute predictors in model
   for(int i = 0; i < data.size(); i++) {
      float obs = data[i].first;
      std::vector<float> ens = data[i].second;
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
      ss << "CalibratorRegression: Cannot train regression, no valid data";
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
}

std::string CalibratorRegression::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c regression", "Applies polynomial regression equation to forecasts: newForecast = a + b * forecast + c * forecast^2 ... . A parameter file is required with the values [a b c ... ]") << std::endl;
   ss << Util::formatDescription("   order=1", "What order is the regression? 0th order: a; 1st order: a + bx; 2nd order a + bc + cx^2; etc.") << std::endl;
   ss << Util::formatDescription("   intercept=1", "Should the regression include an intercept term? If 1, then yes. Only applied when training.") << std::endl;
   return ss.str();
}
