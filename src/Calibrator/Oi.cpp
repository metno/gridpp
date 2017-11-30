#include "Oi.h"
#include "../Util.h"
#include "../Parameters.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <math.h>

CalibratorOi::CalibratorOi(Variable iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mVerticalDecorrelationLength(100),
      mHorizontalDecorrelationLength(30000),
      mMu(0.25),
      mBiasVariable(""),
      mMaxLocations(20),
      mGamma(0.9) {
   iOptions.getValue("bias", mBiasVariable);
   iOptions.getValue("d", mHorizontalDecorrelationLength);
   iOptions.getValue("maxLocations", mMaxLocations);
}

bool CalibratorOi::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   // Check if this method can be applied
   bool hasValidGridpoint = false;
   for(int y = 0; y < nY; y++) {
      for(int x = 0; x < nX; x++) {
         if(Util::isValid(lats[y][x]) && Util::isValid(lons[y][x]) && Util::isValid(elevs[y][x])) {
            hasValidGridpoint = true;
         }
      }
   }
   if(!hasValidGridpoint) {
      Util::warning("There are no gridpoints with valid lat/lon/elev values. Skipping oi...");
      return false;
   }

   std::vector<Location> obsLocations = iParameterFile->getLocations();
   if(iParameterFile->getNumParameters() != 2) {
      std::stringstream ss;
      ss << "Parameter file has " << iParameterFile->getNumParameters() << " parameters, not 2";
      Util::error(ss.str());
   }

   // Create nearest neighbour tree
   std::vector<float> ci(obsLocations.size());
   std::vector<float> obs(obsLocations.size());
   for(int i = 0; i < obsLocations.size(); i++) {
      Parameters parameters = iParameterFile->getParameters(0, obsLocations[i]);
      obs[i] = parameters[0];
      ci[i] = parameters[1];
   }
   float sigma = 2;

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      FieldPtr bias;
      if(iFile.hasVariable(mBiasVariable))
         bias = iFile.getField(mBiasVariable, t);
      FieldPtr accum = iFile.getEmptyField(0);
      FieldPtr weights = iFile.getEmptyField(0);

      #pragma omp parallel for
      for(int y = 0; y < nY; y++) {
         std::cout << y << std::endl;
         double time_s = Util::clock();
         float totalLocations = 0;
         for(int x = 0; x < nX; x++) {
            std::stringstream ss;
            const Location currLocation(lats[y][x], lons[y][x]);
            std::vector<int> useLocations;
            useLocations.reserve(obsLocations.size());

            // Compute observation error covariance matrix

            for(int i = 0; i < obsLocations.size(); i++) {
               float dist = Util::getDistance(obsLocations[i].lat(), obsLocations[i].lon(), lats[y][x], lons[y][x], true);
               if(useLocations.size() > 20)
                  break;
               if(dist < mHorizontalDecorrelationLength) {
                  useLocations.push_back(i);
               }
            }
            totalLocations += useLocations.size()*useLocations.size();

            /*
            int S = useLocations.size();
            ss << "Computing " << y << " " << x << " " << S;
            Util::info(ss.str());
            vec2 Rmatrix(S);
            // vec2 Rmatrix_inv(S);
            for(int i = 0; i < S; i++) {
               Rmatrix[i].resize(S, 0);
               Rmatrix[i][i] = sigma * ci[i];
               //Rmatrix_inv[i].resize(S, 0);
               //Rmatrix_inv[i][i] = 1 / (sigma * ci[i]);
            }
            vec2 Rmatrix_inv = Util::inverse(Rmatrix);

            // Location nearestLocation(Util::MV, Util::MV);
            // iParameterFile->getNearestLocation(t, currLocation, nearestLocation);
            const Parameters parameters = iParameterFile->getParameters(t, currLocation);
            for(int e = 0; e < nEns; e++) {
               (*field)(y, x, e) = parameters[0];
            }
  */
         }
         double time_e = Util::clock();
         std::cout << time_e - time_s << " " << totalLocations / nX << " " << (time_e - time_s) / (totalLocations / nX)<< std::endl;
         abort();
      }
   }
   return true;
}

float CalibratorOi::calcCovar(const Location& loc1, const Location& loc2) const {
   float weight = Util::MV;
   bool isValidLoc1 = Util::isValid(loc1.lat()) && Util::isValid(loc1.lon()) && Util::isValid(loc1.elev());
   bool isValidLoc2 = Util::isValid(loc2.lat()) && Util::isValid(loc2.lon()) && Util::isValid(loc2.elev());
   float useApproxDistance = true;
   if(isValidLoc1 && isValidLoc2) {
      float horizDist = Util::getDistance(loc1.lat(), loc1.lon(), loc2.lat(), loc2.lon(), useApproxDistance);
      float vertDist  = fabs(loc1.elev() - loc2.elev());
      if(horizDist >= mHorizontalDecorrelationLength || (Util::isValid(mVerticalDecorrelationLength) && vertDist >= mVerticalDecorrelationLength)) {
         // Outside radius of influence
         return 0;
      }
      float horizWeight = (mHorizontalDecorrelationLength*mHorizontalDecorrelationLength - horizDist*horizDist) /
                          (mHorizontalDecorrelationLength*mHorizontalDecorrelationLength + horizDist*horizDist);
      float vertWeight = 1;
      if(Util::isValid(mVerticalDecorrelationLength)) {
         vertWeight  = (mVerticalDecorrelationLength*mVerticalDecorrelationLength - vertDist*vertDist) /
                             (mVerticalDecorrelationLength*mVerticalDecorrelationLength + vertDist*vertDist);
      }
      weight = horizWeight * vertWeight;
   }
   return weight;
}

std::string CalibratorOi::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c oi","Spreads bias in space by using kriging. A parameter file is required, which must have one column with the bias.")<< std::endl;
   return ss.str();
}
