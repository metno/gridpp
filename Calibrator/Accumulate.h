#ifndef CALIBRATOR_ACCUMULATE_H
#define CALIBRATOR_ACCUMULATE_H
#include "Calibrator.h"

// Accumlates a certain variable
class CalibratorAccumulate : public Calibrator {
   public:
      CalibratorAccumulate(Variable::Type iFrom, Variable::Type iTo);
   private:
      void calibrateCore(File& iFile) const;
      Variable::Type mFrom;
      Variable::Type mTo;
};
#endif
