#ifndef CALIBRATOR_MASK_H
#define CALIBRATOR_MASK_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Sets gridpoints that are far away from parameter points to missing
class CalibratorMask : public Calibrator {
   public:
      CalibratorMask(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "mask";};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      Variable::Type mVariable;
      bool mUseNearestOnly;
      bool mKeep;
};
#endif
