#ifndef CALIBRATOR_DIAGNOSE_DEWPOINT_H
#define CALIBRATOR_DIAGNOSE_DEWPOINT_H
#include "Calibrator.h"
#include "../Variable.h"

class CalibratorDiagnoseDewpoint : public Calibrator {
   public:
      CalibratorDiagnoseDewpoint(const Variable& iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "diagnoseDrewpoint";};
      bool requiresParameterFile() const { return false;};
      // Temperature in K, RH in [0,1], Returns TD in K
      static float RH2dewpoint(float iTemperature, float iRelativeHumidity);
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      std::string mTemperature;
      std::string mRH;
};
#endif
