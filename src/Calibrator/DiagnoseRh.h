#ifndef CALIBRATOR_DIAGNOSE_RH_H
#define CALIBRATOR_DIAGNOSE_RH_H
#include "Calibrator.h"
#include "../Variable.h"

class CalibratorDiagnoseRh : public Calibrator {
   public:
      CalibratorDiagnoseRh(const Variable& iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "diagnoseDrewpoint";};
      bool requiresParameterFile() const { return false;};
      // Temperatures in K, Returns RH in [0,1]
      static float dewpoint2RH(float iTemperature, float iDewPointTemperature);
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      static float mEwt[41];
      std::string mTemperature;
      std::string mDewpoint;
};
#endif
