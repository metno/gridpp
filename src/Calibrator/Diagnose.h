#ifndef CALIBRATOR_DIAGNOSE_H
#define CALIBRATOR_DIAGNOSE_H
#include "Calibrator.h"
#include "../Variable.h"

class CalibratorDiagnose : public Calibrator {
   public:
      CalibratorDiagnose(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const override {return "diagnose";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
      Variable::Type mOutputVariable;
      // std::vector<Variable::Type> mDiagVariables;
};
#endif
