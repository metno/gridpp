#ifndef CALIBRATOR_NEIGHBOURHOOD_H
#define CALIBRATOR_NEIGHBOURHOOD_H
#include "Calibrator.h"
#include "../Variable.h"
#include "../Util.h"
#include "gridpp.h"

class ParameterFile;
class Parameters;

//! Applies a statistical operator to a neighbourhood
class CalibratorNeighbourhood : public Calibrator {
   public:
      CalibratorNeighbourhood(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "neighbourhood";};
      bool requiresParameterFile() const { return false;};
      void calibrateField(const Field& iInput, Field& iOutput, const Parameters* iParameters=NULL) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      int mRadius;
      gridpp::Statistic mStatType;
      std::string mStatTypeName;
      float mQuantile;
      bool mFast;
      int mNumThresholds;
};
#endif
