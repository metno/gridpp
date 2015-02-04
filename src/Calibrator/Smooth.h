#ifndef CALIBRATOR_SMOOTH_H
#define CALIBRATOR_SMOOTH_H
#include "Calibrator.h"
#include "../Variable.h"

class ParameterFile;
class Parameters;

class CalibratorSmooth : public Calibrator {
   public:
      CalibratorSmooth(Variable::Type iVariable);
      static std::string description();
      std::string name() const {return "smooth";};
      void setSmoothRadius(int iNumPoints);
      int  getSmoothRadius() const;
   private:
      bool calibrateCore(File& iFile) const;
      Variable::Type mVariable;
      int mSmoothRadius;
};
#endif
