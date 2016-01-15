#ifndef CALIBRATOR_CLOUD_H
#define CALIBRATOR_CLOUD_H
#include "Calibrator.h"
#include "../Variable.h"

//! Ensures that if a member has precip it also has full cloud cover
class CalibratorCloud : public Calibrator {
   public:
      CalibratorCloud(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "cloud";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      Variable::Type mPrecipType;
      Variable::Type mCloudType;
};
#endif
