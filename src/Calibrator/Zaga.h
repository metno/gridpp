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
      CalibratorZaga(const ParameterFile* iParameterFile, Variable::Type iMainPredictor);
      //! Get probability mass at 0 mm (i.e probability of no precipitation)
      //! If any input has missing values, the end result is missing
      static float getP0(float iEnsMean, float iEnsFrac, Parameters& iParameters);
      //! Get Precipitation amount corresponding to quantile
      //! If any input has missing values, the end result is missing
      static float getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, Parameters& iParameters);
      //! Set threshold between precip/no-precip used when computing what fraction of
      //! members have precip.
      //! @param must be valid and >= 0
      void setFracThreshold(float iFraction);
      float getFracThreshold() {return mFracThreshold;};

      static std::string description();
      std::string name() const {return "zaga";};
   private:
      const ParameterFile* mParameterFile;
      bool calibrateCore(File& iFile) const;
      static const float mMaxEnsMean;
      //! What precip threshold should be used to count members with no precip?
      float mFracThreshold;
      Variable::Type mMainPredictor;
};
#endif
