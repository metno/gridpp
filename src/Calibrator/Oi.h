#ifndef CALIBRATOR_OI_H
#define CALIBRATOR_OI_H
#include <armadillo>
#include "Calibrator.h"
#include "../ParameterFile/ParameterFile.h"
class Obs;
class Forecast;
class Parameters;

class CalibratorOi : public Calibrator {
   public:
      CalibratorOi(Variable iVariable, const Options& iOptions);
      float calcRho(float iHdist, float iVdist, float iLdist) const;
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
      std::string mDeltaVariable;
      std::string mNumVariable;
      float mC;
      float mEpsilon;
      bool mObsOnly;
      int mMinObs;
      int mX;
      int mY;
      float mMinRho;
      float mMaxBytes;
      bool mUseMeanBias;
      int mMinValidEns;
      bool mUseRho;
      bool mUseBug;
      int mMethod;
      bool mSaveDiff;
      float mElevGradient;
      int mTest;
      float mWMin;
      float mNewDeltaVar;
      float mMaxElevDiff;
      bool mDiagnose;
      bool mUseEns;
      typedef arma::mat mattype;
      typedef arma::vec vectype;
      typedef arma::cx_mat cxtype;
      float calcDelta(float iOldDelta, const vec2& iY) const;
};
#endif
