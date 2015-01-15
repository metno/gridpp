#ifndef CALIBRATOR_ACCUMULATE_H
#define CALIBRATOR_ACCUMULATE_H
#include "Calibrator.h"

// Accumlates a certain variable
class CalibratorAccumulate : public Calibrator {
   public:
      CalibratorAccumulate(const ParameterFile& iParameterFile, Variable::Type iType);
   private:
      void calibrateCore(File& iFile) const;
      Variable::Type mType;
};
#endif
