#ifndef CALIBRATOR_DIAGNOSE_WIND_SPEED_H
#define CALIBRATOR_DIAGNOSE_WIND_SPEED_H
#include "Calibrator.h"
#include "../Variable.h"

class CalibratorDiagnoseWindSpeed : public Calibrator {
   public:
      CalibratorDiagnoseWindSpeed(const Variable& iVariable, const Options& iOptions);
      std::string name() const {return "diagnoseWindSpeed";};
      bool requiresParameterFile() const { return false;};
      static std::string description();
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      std::string mXVariable;
      std::string mYVariable;
};
#endif
