#include "Calibration.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "Util.h"

Calibration::Calibration(const ParameterFile& iParameterFile):
      mParameterFile(iParameterFile) {
}

void CalibrationPrecip::calibrate(const DataFile& iInput, DataFile& iOutput) const {
   int nLat = iInput.getNumLat();
   int nLon = iInput.getNumLon();
   int nEns = iInput.getNumEns();
   int nTime = iInput.getNumTime();
   std::vector<Field*> precipCals(nTime);
   std::vector<Field*> cloudCals(nTime);
   for(int t = 0; t < nTime; t++) {
      precipCals[t] = &iOutput.getEmptyField();
      cloudCals[t] = &iOutput.getEmptyField();
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Parameters parameters = mParameterFile.getParameters(t);
      const Field& precip = iInput.getField(Variable::Precip, t);
      const Field& cloud  = iInput.getField(Variable::Cloud, t);

      // Initialize
      Field& precipCal = *precipCals[t];
      Field& cloudCal  = *cloudCals[t];

      // Parallelizable
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            // Compute model variables
            float ensMean = 0;
            float ensFrac = 0;
            int counter = 0;
            for(int e = 0; e < nEns; e++) {
               ensMean += precip[i][j][e];
               ensFrac += (precip[i][j][e] == 0);//<= mFracThreshold);
               counter++;
            }
            ensMean = ensMean / counter;
            ensFrac = ensFrac / counter;

            // Check inputs to model
            if(ensMean > mMaxEnsMean) {
               ensMean = mMaxEnsMean;
            }

            // TODO: Figure out which cloudless members to use
            float p0 = getP0(ensMean, ensFrac, parameters);
            // Calibrate
            std::vector<std::pair<float,int> > pairs(nEns);
            std::vector<float> valuesCal(nEns);
            for(int e = 0; e < nEns; e++) {
               float quantile = ((float) e+0.5)/nEns;
               float valueCal   = getInvCdf(quantile, ensMean, ensFrac, t, parameters);
               float valueUncal = precip[i][j][e];
               valuesCal[e] = valueCal;
               pairs[e].first = valueUncal;
               pairs[e].second = e;
            }
            // TODO: Resort members
            std::sort(pairs.begin(), pairs.end(), sort_pair_first<float,int>());
            for(int e = 0; e < nEns; e++) {
               int ei = pairs[e].second;
               float valueCal = valuesCal[e];
               precipCal[i][j][ei] = valueCal;
            }

            // Fix cloud cover if rain
            for(int e = 0; e < nEns; e++) {
               float precip = precipCal[i][j][e];
               cloudCal[i][j][e]  = cloud[i][j][e];

               if(precip > 0 && cloud[i][j][e] == 0) {
                  cloudCal[i][j][e]  = 1;
               }
            }
         }
      }
   }

   // Accumulate precipitation
   for(int t = 1; t < nTime; t++) {
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               (*precipCals[t])[i][j][e] = (*precipCals[t])[i][j][e] + (*precipCals[t-1])[i][j][e];
            }
         }
      }
   }

   // Add to file
   for(int t = 0; t < nTime; t++) {
      iOutput.addField(*precipCals[t], Variable::PrecipAcc, t);
      iOutput.addField(*cloudCals[t], Variable::Cloud, t);
   }
}

float CalibrationPrecip::getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, int iTime, Parameters& iParameters) {
   if(iEnsMean < 0)
      abort();
   float P0 = getP0(iEnsMean, iEnsFrac, iParameters);
   if(iQuantile < P0)
      return 0;

   float mua = iParameters[0];
   float mub = iParameters[1];
   float muc = iParameters[2];
   float sa = iParameters[3];
   float sb = iParameters[4];
   float sc = iParameters[5];

   float quantileCont = (iQuantile-P0)/(1-P0);
   // Code to get CDF from Gamma distribution
   // Parameters in R:
   float mu    = exp(mua + mub * pow(iEnsMean, 1.0/3) + muc * iTime);
   float sigma = exp(sa + sb * iEnsMean + sc * iTime);

   if(sigma == 0)
      abort();

   // Parameters in boost and wikipedia
   float shape = 1/(sigma*sigma); // k
   float scale = sigma*sigma*mu;  // theta
   // std::cout << mu << " " << sigma << " " << P0 << " " << shape << " " << scale << std::endl;
   boost::math::gamma_distribution<> dist(shape, scale);
   float value = boost::math::quantile(dist, quantileCont);
   return value;
}
float CalibrationPrecip::getP0(float iEnsMean, float iEnsFrac, Parameters& iParameters) {
   float a = iParameters[6];
   float b = iParameters[7];
   float c = iParameters[8];
   float logit = a + b * iEnsMean + c * iEnsFrac;
   float P0 = invLogit(logit);
   return P0;
}
float Calibration::logit(float p) {
   return log(p/(1-p));
}
float Calibration::invLogit(float x) {
   return exp(x)/(exp(x)+1);
}

CalibrationPrecip::CalibrationPrecip(const ParameterFile& iParameterFile):
      Calibration(iParameterFile),
      mFracThreshold(0.5) {
}
