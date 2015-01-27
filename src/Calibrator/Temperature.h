#ifndef CALIBRATOR_TEMPERATURE_H
#define CALIBRATOR_TEMPERATURE_H
#include "Calibrator.h"
#include "../ParameterFile.h"
class CalibratorTemperature : public Calibrator {
   public:
      CalibratorTemperature(const ParameterFile& iParameterFile);
      bool calibrateCore(const File& iInput, File& iOutput) const;

};
#endif
