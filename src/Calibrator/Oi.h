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
      float calcRho(float iHdist, float iVdist) const;
      static std::string description();
      std::string name() const {return "oi";};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      float mVLength;
      float mHLength;
      float mMu;
      float mGamma;
      std::string mBiasVariable;
      int mMaxLocations;
      bool mSort;
      float mSigma;
      float mDelta;
      std::string mNumVariable;
      bool mObsOnly;
      int mMinObs;
      int mX;
      int mY;
      float mMinRho;
      float mMaxBytes;
      bool mUseMeanBias;
};
#endif
