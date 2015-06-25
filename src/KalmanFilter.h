#ifndef KALMANFILTER_H
#define KALMANFILTER_H
#include "Variable.h"
class Obs;
class Forecast;
class Parameters;
class File;

class KalmanFilter {
   public:
      KalmanFilter(Variable::Type iVariable, float iRatio);
      //! Returns the next bias
      float calcBias(const Obs& iObs, const Forecast& iForecast, const Parameters& iParameters) const;
      float calcBias(const Parameters& iParameters) const;
      //! Updates parameters
      Parameters update(const Obs& iObs, const Forecast& iForecast, const Parameters& iOldParameters) const;
      Parameters update(float iObs, float iForecast, const Parameters& iOldParameters) const;

   private:
      Variable::Type mVariable;
      float mRatio;
      static float mVarVarV;
      static float mVarVarW;
      static float mMaxP;
};
#endif
