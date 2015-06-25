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
      ~CalibratorKriging();
      static std::string description();
      float calcWeight(const Location& loc1, const Location& loc2) const;
      enum Type {
         TypeCressman = 10,
         TypeBarnes   = 20
      };
   private:
      bool calibrateCore(File& iFile) const;

      Variable::Type mVariable;
      float mRadius;
      float mMaxElevDiff;
      float mEfoldDist;
      const ParameterFileSpatial* mParameterFile;
      std::string name() const {return "kriging";};
      File* mPrevious;
      Variable::Type mAuxVariable;
      float mLowerThreshold;
      float mUpperThreshold;
      Type mKrigingType;
      bool mUseApproxDistance;
      Util::Operator mOperator;
};
#endif
