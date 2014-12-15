#include "Calibration.h"

Calibration::Calibration(const ParameterFile& iParameterFile):
      mParameterFile(iParameterFile),
      mPrecipName("precipitation"), 
      mCloudName("cloudCover") {
}

void Calibration::calibrate(const DataFile& iInput, DataFile& iOutput) const {
   std::vector<int> offsets;
   int nX = iInput.getNumX();
   int nY = iInput.getNumY();
   int nE = iInput.getNumE();
   // Loop over offsets
   for(int t = 0; t < offsets.size(); t++) {
      Parameters parameters = mParameterFile.getParameters(t);
      const Field& precip = iInput.getField(mPrecipName, t);
      const Field& cloud  = iInput.getField(mCloudName, t);

      // Initialize
      Field precipCal;
      Field cloudCal;

      // Parallelizable
      for(int i = 0; i < nX; i++) {
         for(int j = 0; j < nY; j++) {
            // Compute model variables
            float ensMean = 0;
            float ensFrac = 0;
            int counter = 0;
            for(int e = 0; e < nE; e++) {
               ensMean += precip[i][j][e];
               ensFrac += (precip[i][j][e] > 0);
               counter++;
            }
            ensMean = ensMean / counter;
            ensFrac = ensFrac / counter;

            // Check inputs to model
            if(ensMean > mMaxEnsMean) {
               ensMean = mMaxEnsMean;
            }

            // TODO: Figure out which cloudless members to use
            float p0 = getP0(ensMean, ensMean, parameters);
            // Calibrate
            for(int e = 0; e < nE; e++) {
               // TODO: Resort members
               float quantile = ((float) e+0.5)/nE;
               float value = getInvCdf(quantile, ensMean, ensFrac, parameters);
               precipCal[i][j][e] = value;

               // Fix cloud cover if rain
               if(value > 0 && cloud[i][j][e] == 0) {
                  cloudCal[i][j][e]  = value;
               }
            }
         }
      }
      iOutput.setField(precipCal, mPrecipName, t);
      iOutput.setField(cloudCal, mCloudName, t);
   }
}

float Calibration::getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, Parameters& iParameters) const {
   // Code to get CDF from Gamma distribution
   return 1;
}
float Calibration::getP0(float iEnsMean, float iEnsFrac, Parameters& iParameters) const {
   return 0.3;
}
