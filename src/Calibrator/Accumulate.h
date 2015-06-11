#ifndef CALIBRATOR_ACCUMULATE_H
#define CALIBRATOR_ACCUMULATE_H
#include "Calibrator.h"
#include "../Variable.h"

// Accumlates a certain variable
class CalibratorAccumulate : public Calibrator {
   public:
      CalibratorAccumulate(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "accumulate";};
   private:
      bool calibrateCore(File& iFile) const;
      Variable::Type mInputVariable;
      Variable::Type mOutputVariable;
};
#endif
