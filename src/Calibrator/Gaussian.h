#ifndef CALIBRATOR_GAUSSIAN_H
#define CALIBRATOR_GAUSSIAN_H
#include "Calibrator.h"
#include "../Variable.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

class ParameterFile;
class Parameters;

//! Ensemble calibration using zero-adjusted gamma distribution. Its predictors are:
//! - ensemble mean
//! - ensemble fraction
//! Designed for precip
class CalibratorGaussian : public Calibrator {
   public:
      CalibratorGaussian(Variable::Type iMainPredictor, const Options& iOptions);
      static float getInvCdf(float iQuantile, float iEnsMean, float iEnsSpread, const Parameters& iParameters);
      static float getCdf(float iThreshold, float iEnsMean, float iEnsSpread, const Parameters& iParameters);
      static float getPdf(float iThreshold, float iEnsMean, float iEnsSpread, const Parameters& iParameters);

      int   getNeighbourhoodSize() {return mNeighbourhoodSize;};

      static std::string description();
      std::string name() const override {return "gaussian";};
      Parameters train(const std::vector<ObsEns>& iData) const override;
   private:
      static double my_f(const gsl_vector *v, void *params);
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
      Variable::Type mMainPredictor;
      int  mNeighbourhoodSize;
      float mLogLikelihoodTolerance;
      static const int mNumParameters = 2;
};
#endif
