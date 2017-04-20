#ifndef CALIBRATOR_DIAGNOSE_H
#define CALIBRATOR_DIAGNOSE_H
#include "Calibrator.h"
#include "../Variable.h"

class CalibratorDiagnose : public Calibrator {
   public:
      CalibratorDiagnose(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "diagnose";};
      bool requiresParameterFile() const { return false;};
      // Temperatures in K, Returns RH in [0,1]
      static float dewpoint2RH(float iTemperature, float iDewPointTemperature);
      // Temperature in K, RH in [0,1], Returns TD in K
      static float RH2dewpoint(float iTemperature, float iRelativeHumidity);
   private:
      // Data needed for dewpoint to RH conversion
      static float mEwt[41];
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      Variable::Type mOutputVariable;
      // std::vector<Variable::Type> mDiagVariables;
};
#endif
