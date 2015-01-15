#ifndef CALIBRATOR_CLOUD_H
#define CALIBRATOR_CLOUD_H
#include "Calibrator.h"

// Ensure that if a member has precip it also has full cloud cover
class CalibratorCloud : public Calibrator {
   public:
      CalibratorCloud(const ParameterFile& iParameterFile, Variable::Type iPrecip=Variable::Precip, Variable::Type iCloud=Variable::Cloud);
   private:
      void calibrateCore(DataFile& iFile) const;
      Variable::Type mPrecipType;
      Variable::Type mCloudType;
};
#endif
