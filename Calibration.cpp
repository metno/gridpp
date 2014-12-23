#include "Calibration.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "Util.h"

Calibration::Calibration(const ParameterFile& iParameterFile):
      mParameterFile(iParameterFile) {
}

CalibrationPrecip::CalibrationPrecip(const ParameterFile& iParameterFile):
      Calibration(iParameterFile),
      mFracThreshold(0.5) {
}

void CalibrationPrecip::calibrate(const DataFile& iInput, DataFile& iOutput) const {
   int nLat = iInput.getNumLat();
   int nLon = iInput.getNumLon();
   int nEns = iInput.getNumEns();
   int nTime = iInput.getNumTime();

   // Initialize calibrated fields in the output file
   std::vector<Field*> precipCals(nTime);
   std::vector<Field*> cloudCals(nTime);
   for(int t = 0; t < nTime; t++) {
      precipCals[t] = &iOutput.getEmptyField();
      cloudCals[t] = &iOutput.getEmptyField();
   }

   int numInvalidRaw = 0;
   int numInvalidCal = 0;

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Parameters parameters = mParameterFile.getParameters(t);
      const Field& precip = iInput.getField(Variable::Precip, t);
      const Field& cloud  = iInput.getField(Variable::Cloud, t);

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
            bool isValid = true;
            // Check if current ensemble for this gridpoint/time has any missing
            // values. If so, don't calibrate the ensemble.
            for(int e = 0; e < nEns; e++) {
               float value = precip[i][j][e];
               if(!Util::isValid(value)) {
                  isValid = false;
                  break;
               }
               ensMean += value;
               ensFrac += (value <= mFracThreshold);
               counter++;
            }
            if(isValid) {
               ensMean = ensMean / counter;
               ensFrac = ensFrac / counter;

               // Limit the input to the calibration model to prevent it
               // from creating very extreme values.
               if(ensMean > mMaxEnsMean) {
                  ensMean = mMaxEnsMean;
               }

               // TODO: Figure out which cloudless members to use. Ideally, if more members
               // need precip, we should pick members that already have clouds, so that we minimize
               // our effect on the cloud cover field.

               // Calibrate
               std::vector<std::pair<float,int> > pairs(nEns);
               std::vector<float> valuesCal(nEns);
               bool isValid = true;
               for(int e = 0; e < nEns; e++) {
                  float quantile = ((float) e+0.5)/nEns;
                  float valueCal   = getInvCdf(quantile, ensMean, ensFrac, t, parameters);
                  float valueUncal = precip[i][j][e];
                  if(!Util::isValid(valueCal)) {
                     isValid = false;
                     break;
                  }
                  valuesCal[e] = valueCal;
                  pairs[e].first = valueUncal;
                  pairs[e].second = e;
               }
               if(isValid) {
                  // Sort values so that the rank of a member is the same before and after calibration
                  std::sort(pairs.begin(), pairs.end(), sort_pair_first<float,int>());
                  for(int e = 0; e < nEns; e++) {
                     int ei = pairs[e].second;
                     float valueCal = valuesCal[e];
                     precipCal[i][j][ei] = valueCal;
                  }

                  // Turn on clouds if needed, i.e don't allow a member to
                  // have precip without cloud cover.
                  for(int e = 0; e < nEns; e++) {
                     float precip = precipCal[i][j][e];
                     cloudCal[i][j][e]  = cloud[i][j][e];

                     if(precip > 0 && cloud[i][j][e] == 0) {
                        cloudCal[i][j][e]  = 1;
                     }
                  }
               }
               else {
                  numInvalidCal++;
                  // Calibration produced some invalid members. Revert to the raw values.
                  for(int e = 0; e < nEns; e++) {
                     precipCal[i][j][e] = precip[i][j][e];
                  }
               }
            }
            else {
               numInvalidRaw++;
               // One or more members are missing, don't calibrate
               for(int e = 0; e < nEns; e++) {
                  precipCal[i][j][e] = precip[i][j][e];
               }
            }
         }
      }
   }
   if(numInvalidRaw > 0) {
      std::stringstream ss;
      ss << "File '" << iInput.getFilename() << "' has " << numInvalidRaw
         << " missing ensembles, out of " << nTime * nLat * nLon << ".";
      Util::warning(ss.str());
   }
   if(numInvalidCal > 0) {
      std::stringstream ss;
      ss << "Calibration produced " << numInvalidCal
         << " invalid ensembles, out of " << nTime * nLat * nLon << ".";
      Util::warning(ss.str());
   }

   // Accumulate precipitation
   for(int t = 1; t < nTime; t++) {
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               float previous = (*precipCals[t-1])[i][j][e];
               float current = (*precipCals[t])[i][j][e];
               if(Util::isValid(current) && Util::isValid(previous)) {
                  (*precipCals[t])[i][j][e] = current + previous;
               }
               else
                  (*precipCals[t])[i][j][e] = Util::MV;
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
   if(iQuantile == 0)
      return 0;

   if(iQuantile < 0 || iQuantile >= 1) {
      Util::warning("Quantile must be in the interval [0,1)");
      return Util::MV;
   }
   if(!Util::isValid(iEnsMean) || !Util::isValid(iEnsFrac))
      return Util::MV;

   if(iEnsMean < 0 || iEnsFrac < 0 || iEnsFrac > 1 || iTime < 0)
      return Util::MV;

   // Check if we are in the discrete mass
   float P0 = getP0(iEnsMean, iEnsFrac, iParameters);
   if(!Util::isValid(P0))
      return Util::MV;
   if(iQuantile < P0)
      return 0;

   float mua = iParameters[0];
   float mub = iParameters[1];
   float muc = iParameters[2];
   float sa  = iParameters[3];
   float sb  = iParameters[4];
   float sc  = iParameters[5];

   float quantileCont = (iQuantile-P0)/(1-P0);
   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = exp(mua + mub * pow(iEnsMean, 1.0/3) + muc * iTime);
   float sigma = exp(sa + sb * iEnsMean + sc * iTime);

   if(mu <= 0 || sigma <= 0)
      return Util::MV;
   if(!Util::isValid(mu) || !Util::isValid(sigma))
      return Util::MV;

   // Parameters in boost and wikipedia
   float shape = 1/(sigma*sigma); // k
   float scale = sigma*sigma*mu;  // theta
   if(!Util::isValid(scale) || !Util::isValid(shape))
      return Util::MV;

   // std::cout << mu << " " << sigma << " " << P0 << " " << shape << " " << scale << std::endl;
   boost::math::gamma_distribution<> dist(shape, scale);
   float value = boost::math::quantile(dist, quantileCont);
   if(!Util::isValid(value))
      return Util::MV;
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
void CalibrationPrecip::setFracThreshold(float iFraction) {
   if(iFraction < 0 || iFraction > 1) {
      std::stringstream ss;
      ss << "CalibrationPrecip: fraction threshold (" << iFraction << ") must be between 0 and 1, inclusive.";
      Util::error(ss.str());
   }
   mFracThreshold = iFraction;
}
float Calibration::logit(float p) {
   return log(p/(1-p));
}
float Calibration::invLogit(float x) {
   return exp(x)/(exp(x)+1);
}

