#include "KalmanFilter.h"
#include "Util.h"
#include "Parameters.h"
#include "File/File.h"
#include <math.h>
float KalmanFilter::mVarVarV = 1;
float KalmanFilter::mVarVarW = 0.0005;
float KalmanFilter::mMaxP    = 10000;

KalmanFilter::KalmanFilter(Variable::Type iVariable, float iRatio) :
      mVariable(iVariable),
      mRatio(iRatio) {
}

float KalmanFilter::calcBias(const Obs& iObs, const Forecast& iForecast, const Parameters& iParameters) const {
   float kalmanGain              = iParameters[4];
   float yesterdaysError         = iParameters[5];
   float yesterdaysBiasEstimate  = iParameters[7];
   float todaysBiasEstimate = 0;

   if(Util::isValid(yesterdaysError))
      todaysBiasEstimate = yesterdaysBiasEstimate + kalmanGain*(yesterdaysError - yesterdaysBiasEstimate);
   else
      todaysBiasEstimate = yesterdaysBiasEstimate;

   return todaysBiasEstimate;
}
float KalmanFilter::calcBias(const Parameters& iParameters) const {
   float kalmanGain              = iParameters[4];
   float yesterdaysError         = iParameters[5];
   float yesterdaysBiasEstimate  = iParameters[7];
   float todaysBiasEstimate = 0;

   if(Util::isValid(yesterdaysError))
      todaysBiasEstimate = yesterdaysBiasEstimate + kalmanGain*(yesterdaysError - yesterdaysBiasEstimate);
   else
      todaysBiasEstimate = yesterdaysBiasEstimate;

   return todaysBiasEstimate;
}
Parameters KalmanFilter::update(const Obs& iObs, const Forecast& iForecast, const Parameters& iOldParameters) const {
   std::vector<float> values;
   
   return Parameters(values);
}

Parameters KalmanFilter::update(float iObs, float iForecast, const Parameters& iOldParameters) const {
   std::vector<float> parameters(8, Util::MV);
   if(iOldParameters.size() != 8) {
      parameters[0] = 1000;
      parameters[1] = 0;
      parameters[2] = 1;
      parameters[3] = 1;
      parameters[4] = 0;
      parameters[5] = 0;
      parameters[6] = 0;
      parameters[7] = 0;
   }
   else {
      float pVarV             = iOldParameters[0];
      float kalmanGainVar     = iOldParameters[1];
      float varV              = iOldParameters[2];
      float p                 = iOldParameters[3];
      float kalmanGain        = iOldParameters[4];
      float lastError         = iOldParameters[5];
      // float previousLastError = iOldParameters[6];
      float biasEstimate      = iOldParameters[7];

      assert(pVarV > 0);
      assert(varV > 0);
      assert(p > 0);
      assert(kalmanGainVar >= 0 && kalmanGainVar <= 1);
      assert(kalmanGain    >= 0 && kalmanGain    <= 1);

      // Update
      float error = Util::MV;
      if(Util::isValid(iObs) && Util::isValid(iForecast)) {
         error = iForecast - iObs;
         // Set currPvarV currKalmanGainVar, currVarV
         if(Util::isValid(lastError)) {
            assert(pVarV + mVarVarW + mVarVarV > 0);

            kalmanGainVar = (pVarV + mVarVarW) / (pVarV + mVarVarW + mVarVarV);
            pVarV         = (pVarV + mVarVarW) * (1 - kalmanGainVar);

            varV = varV + kalmanGainVar*(pow((error - lastError),2)/(2 + mRatio) - varV);
         }
         float varW = mRatio * varV;

         kalmanGain = (p + varW) / (p + varW + varV);
         p = (p + varW)*(1 - kalmanGain);
         assert(p + varW + varV > 0);

         biasEstimate = biasEstimate + kalmanGain*(error - biasEstimate);
      }
      else {
         error = Util::MV;
         p = p + mRatio * varV;

      }
      if(p > mMaxP) {
         p = mMaxP;
      }

      parameters[0] = pVarV;
      parameters[1] = kalmanGainVar;
      parameters[2] = varV;
      parameters[3] = p;
      parameters[4] = kalmanGain;
      parameters[5] = error;
      parameters[6] = lastError;
      parameters[7] = biasEstimate;
   }

   return Parameters(parameters);
}
