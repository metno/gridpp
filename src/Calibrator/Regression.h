#ifndef CALIBRATOR_REGRESSION_H
#define CALIBRATOR_REGRESSION_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;
class TrainingData;

//! Applies polynomial regression to forecasts
class CalibratorRegression : public Calibrator {
   public:
      CalibratorRegression(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "regression";};
      Parameters train(const std::vector<ObsEns>& iData) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      Variable::Type mVariable;
      int mOrder;
      bool mIntercept;
};
#endif
