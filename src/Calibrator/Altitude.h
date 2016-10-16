#ifndef CALIBRATOR_ALTITUDE_H
#define CALIBRATOR_ALTITUDE_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Changes the altitudes in iFile to the altitudes in the parameter file
class CalibratorAltitude : public Calibrator {
   public:
      CalibratorAltitude(const Options& iOptions);
      static std::string description();
      std::string name() const {return "altitude";};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
};
#endif
