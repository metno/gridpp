#ifndef CALIBRATOR_REGRESSION_H
#define CALIBRATOR_REGRESSION_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Applies polynomial regression to forecasts
class CalibratorRegression : public Calibrator {
   public:
      CalibratorRegression(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const override {return "regression";};
      Parameters train(const std::vector<ObsEns>& iData) const override;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
      Variable::Type mVariable;
      int mOrder;
      bool mIntercept;
};
#endif
