#ifndef CALIBRATOR_SMOOTH_H
#define CALIBRATOR_SMOOTH_H
#include "Calibrator.h"
#include "../Variable.h"

class File;
class ParameterFile;
class Parameters;

class CalibratorSmooth : public Calibrator {
   public:
      CalibratorSmooth(Variable::Type iMainPredictor);
      static std::string description();
      std::string name() const {return "Smooth";};
      void setSmoothRadius(int iNumPoints);
   private:
      void calibrateCore(File& iFile) const;
      Variable::Type mMainPredictor;
      int mSmoothRadius;
};
#endif
