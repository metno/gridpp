#ifndef CALIBRATOR_WIND_H
#define CALIBRATOR_WIND_H
#include "Calibrator.h"
#include "../ParameterFile.h"

class CalibratorWind : public Calibrator {
   public:
      CalibratorWind(const ParameterFileRegion& iParameterFileRegion);
   private:
      void calibrateCore(File& iFile) const override;
      static float getDir(float iU, float iV);
      static float getSpeed(float iU, float iV);
      float calibrateLocal(float iSpeed, float iDir, const Parameters& iParameters) const;
      ParameterFileRegion mParameterFile;
};
#endif
