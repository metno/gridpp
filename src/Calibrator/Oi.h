#ifndef CALIBRATOR_OI_H
#define CALIBRATOR_OI_H
#include "Calibrator.h"
#include "../ParameterFile/ParameterFile.h"
class Obs;
class Forecast;
class Parameters;

class CalibratorOi : public Calibrator {
   public:
      CalibratorOi(Variable iVariable, const Options& iOptions);
      static std::string description();
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      float calcCovar(const Location& loc1, const Location& loc2) const;
      std::string name() const {return "oi";};
      float mVerticalDecorrelationLength;
      float mHorizontalDecorrelationLength;
      float mMu;
      float mGamma;
      std::string mBiasVariable;
      int mMaxLocations;
};
#endif
