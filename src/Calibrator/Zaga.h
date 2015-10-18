#ifndef CALIBRATOR_ZAGA_H
#define CALIBRATOR_ZAGA_H
#include "Calibrator.h"
#include "../Variable.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

class ParameterFile;
class Parameters;
class TrainingData;

//! Ensemble calibration using zero-adjusted gamma distribution. Its predictors are:
//! - ensemble mean
//! - ensemble fraction
//! Designed for precip
class CalibratorZaga : public Calibrator {
   public:
      CalibratorZaga(Variable::Type iMainPredictor, const Options& iOptions);
      //! Get probability mass at 0 mm (i.e probability of no precipitation)
      //! If any input has missing values, the end result is missing
      static float getP0(float iEnsMean, float iEnsFrac, const Parameters& iParameters);
      //! Get Precipitation amount corresponding to quantile
      //! If any input has missing values, the end result is missing
      static float getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, const Parameters& iParameters);
      //! Get Precipitation amount corresponding to quantile. If any input has missing values, or
      //! iEnsMean < 0 or iEnsFrac is not in [0,1], a missing value is returned.
      static float getCdf(float iThreshold, float iEnsMean, float iEnsFrac, const Parameters& iParameters);
      static float getPdf(float iThreshold, float iEnsMean, float iEnsFrac, const Parameters& iParameters);

      float getFracThreshold() {return mFracThreshold;};
      bool  getOutputPop() {return mOutputPop;};
      float getPopThreshold() {return mPopThreshold;};
      int   getNeighbourhoodSize() {return mNeighbourhoodSize;};
      float getMaxEnsMean() {return mMaxEnsMean;};

      static std::string description();
      std::string name() const {return "zaga";};
      //! Create parameters based on training data
      Parameters train(const TrainingData& iData, int iOffset) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      static float logLikelihood(float obs, float iEnsMean, float iEnsFrac, const Parameters& iParameters);
      static double my_f(const gsl_vector *v, void *params);
      // static void my_df(const gsl_vector *v, void *params, gsl_vector *df);
      // static void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);
      //! What precip threshold should be used to count members with no precip?
      float mFracThreshold;
      Variable::Type mMainPredictor;
      bool mOutputPop;
      int  mNeighbourhoodSize;
      float mPopThreshold;
      float mPrecipLowQuantile;
      float mPrecipMiddleQuantile;
      float mPrecipHighQuantile;
      float mMaxEnsMean;
      bool m6h;
      float mLogLikelihoodTolerance;
};
#endif
