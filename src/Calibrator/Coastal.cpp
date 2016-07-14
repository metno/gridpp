#include "Coastal.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
#include "../Grid.h"
CalibratorCoastal::CalibratorCoastal(Variable::Type iVariable, const Options& iOptions) :
      Calibrator(iOptions),
      mVariable(iVariable),
      mSearchRadius(3),
      mMinLsmDiff(0.4) {
   iOptions.getValue("searchRadius", mSearchRadius);
   iOptions.getValue("minLsmDiff", mMinLsmDiff);
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
               float minLsm = 2;
               float maxLsm = -1;
               float nearestNeighbour = (*field)(i,j,e);
               // Find the range within a neighbourhood
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

               float a = parameters[0];
               float b = parameters[1];
               float c = parameters[2];
               // If there is not sufficient variability in LSM in the neighbourhood, then
               // force don't use the range in the regression
               float range = 0;
               if(Util::isValid(lowerValue) && Util::isValid(upperValue) && (maxLsm - minLsm > mMinLsmDiff))
                  range = upperValue - lowerValue;

               float value = a + b * nearestNeighbour + c * range;
               (*field)(i,j,e) = value;
            }
         }
      }
   }
   return true;
}

Parameters CalibratorCoastal::train(const std::vector<ObsEnsField>& iData, const Grid& iObsGrid, const Grid& iEnsGrid, int iIobs, int iJobs, int iIens, int iJens) const {
   std::vector<float> y;
   std::vector<std::vector<float> > x(3);
   int nLat = iData[0].second->getNumLat();
   int nLon = iData[0].second->getNumLon();
   vec2 lsm = iEnsGrid.landFractions();

   // Find nearest land and sea point
   float lowerValue = Util::MV;
   float upperValue = Util::MV;
   int Imin = Util::MV;
   int Imax = Util::MV;
   int Jmin = Util::MV;
   int Jmax = Util::MV;
   float minLsm = 2;
   float maxLsm = -1;
   // Find the range within a neighbourhood
   for(int ii = std::max(0, iIens-mSearchRadius); ii <= std::min(nLat-1, iIens+mSearchRadius); ii++) {
      for(int jj = std::max(0, iJens-mSearchRadius); jj <= std::min(nLon-1, iJens+mSearchRadius); jj++) {
         float currLsm = lsm[ii][jj];
         if(currLsm < minLsm) {
            Imin = ii;
            Jmin = jj;
            minLsm = currLsm;
         }
         if(currLsm > minLsm) {
            Imax = ii;
            Jmax = jj;
            maxLsm = currLsm;
         }
      }
   }

   assert(Util::isValid(Imax));
   assert(Util::isValid(Jmax));

   // Compute predictors in model
   for(int i = 0; i < iData.size(); i++) {
      float obs = (*iData[i].first)(iIobs, iJobs, 0);
      std::vector<float> ensMin = (*iData[i].second)(Imin, Jmin);
      float valueMin = Util::calculateStat(ensMin, Util::StatTypeMean);
      std::vector<float> ensMax = (*iData[i].second)(Imax, Jmax);
      float valueMax = Util::calculateStat(ensMax, Util::StatTypeMean);
      std::vector<float> ensNn = (*iData[i].second)(iIens, iJens);
      float valueNn = Util::calculateStat(ensNn, Util::StatTypeMean);
      if(Util::isValid(obs) && Util::isValid(valueNn)) {
         float range = 0;
         if(Util::isValid(valueMin) && Util::isValid(valueMax))
            range = valueMax - valueMin;
         y.push_back(obs);
         x[0].push_back(1);
         x[1].push_back(valueNn);
         x[2].push_back(range);
      }
   }

   if(y.size() <= 0) {
      std::stringstream ss;
      ss << "CalibratorCoastal: Cannot train coastal, no valid data";
      Util::error(ss.str());
   }

   std::vector<float> values = Util::regression(y, x, false);
   Parameters par(values);

   return par;
}

std::string CalibratorCoastal::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c coastal", "Weights land and sea forecasts") << std::endl;
   ss << Util::formatDescription("   searchRadius=3", "") << std::endl;
   ss << Util::formatDescription("   minLsmDiff=0.4", "Minimum difference in LSM in order to use the  land and sea points") << std::endl;
   return ss.str();
}
