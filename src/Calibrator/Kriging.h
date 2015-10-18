#ifndef CALIBRATOR_KRIGING_H
#define CALIBRATOR_KRIGING_H
#include "Calibrator.h"
#include "../ParameterFile/ParameterFile.h"
class Obs;
class Forecast;
class Parameters;

class CalibratorKriging : public Calibrator {
   public:
      CalibratorKriging(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      float calcCovar(const Location& loc1, const Location& loc2) const;
      enum Type {
         TypeCressman = 10,
         TypeBarnes   = 20
      };
      //! Compute the bias at the training point
      Parameters train(const TrainingData& iData, int iOffset) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;

      Variable::Type mVariable;
      float mRadius;
      float mMaxElevDiff;
      float mEfoldDist;
      std::string name() const {return "kriging";};
      File* mPrevious;
      Variable::Type mAuxVariable;
      float mLowerThreshold;
      float mUpperThreshold;
      Type mKrigingType;
      bool mUseApproxDistance;
      Util::Operator mOperator;
      bool mCrossValidate;
};
#endif
