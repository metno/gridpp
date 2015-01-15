#ifndef CALIBRATOR_ZAGA_H
#define CALIBRATOR_ZAGA_H
#include "Calibrator.h"
#include "../Variable.h"

class File;
class ParameterFile;
class Parameters;

// Zero-adjusted gamma distribution, using predictors:
// - ensemble mean
// - ensemble fraction
// Designed for precip
class CalibratorZaga : public Calibrator {
   public:
      CalibratorZaga(const ParameterFile& iParameterFile, Variable::Type iMainPredictor);
      //! Get probability mass at 0 mm (i.e probability of no precipitation)
      static float getP0(float iEnsMean, float iEnsFrac, Parameters& iParameters);
      //! Get Precipitation amount corresponding to quantile
      static float getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, int iTime, Parameters& iParameters);
      //! Set threshold between precip/no-precip used when computing what fraction of
      //! members have precip.
      void setFracThreshold(float iFraction);
   private:
      const ParameterFile& mParameterFile;
      void calibrateCore(File& iFile) const;
      static const float mMaxEnsMean = 100;
      //! What precip threshold should be used to count members with no precip?
      float mFracThreshold;
      int mNeighbourhoodSize;
      Variable::Type mMainPredictor;
};
#endif
