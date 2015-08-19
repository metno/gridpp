#ifndef CALIBRATOR_QNH_H
#define CALIBRATOR_QNH_H
#include "Calibrator.h"
#include "../Variable.h"

//! Creates the QNH field by using surface pressure
class CalibratorQnh : public Calibrator {
   public:
      CalibratorQnh(const Options& iOptions);
      static std::string description();
      std::string name() const {return "qnh";};
      static float calcQnh(float iElev, float iPressure);
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
};
#endif
