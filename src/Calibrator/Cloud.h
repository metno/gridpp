#ifndef CALIBRATOR_CLOUD_H
#define CALIBRATOR_CLOUD_H
#include "Calibrator.h"
#include "../Variable.h"

//! Ensures that if a member has precip it also has full cloud cover
class CalibratorCloud : public Calibrator {
   public:
      CalibratorCloud(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const override {return "cloud";};
      bool requiresParameterFile() const override { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
      Variable::Type mPrecipType;
      Variable::Type mCloudType;
};
#endif
