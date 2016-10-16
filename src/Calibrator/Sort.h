#ifndef CALIBRATOR_SORT_H
#define CALIBRATOR_SORT_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Applies polynomial regression to forecasts
class CalibratorSort : public Calibrator {
   public:
      CalibratorSort(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "sort";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
      Variable::Type mVariable;
};
#endif
