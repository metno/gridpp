#ifndef CALIBRATOR_REGRESSION_H
#define CALIBRATOR_REGRESSION_H
#include "Calibrator.h"
class ParameterFile;

//! Applies polynomial regression to forecasts
class CalibratorRegression : public Calibrator {
   public:
      CalibratorRegression(const ParameterFile* iParameterFile, Variable::Type iVariable);
      static std::string description();
      std::string name() const {return "regression";};
   private:
      bool calibrateCore(File& iFile) const;
      Variable::Type mVariable;
};
#endif
