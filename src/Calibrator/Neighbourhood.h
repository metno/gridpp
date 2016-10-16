#ifndef CALIBRATOR_NEIGHBOURHOOD_H
#define CALIBRATOR_NEIGHBOURHOOD_H
#include "Calibrator.h"
#include "../Variable.h"
#include "../Util.h"

class ParameterFile;
class Parameters;

//! Applies a statistical operator to a neighbourhood
class CalibratorNeighbourhood : public Calibrator {
   public:
      CalibratorNeighbourhood(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const override {return "neighbourhood";};
      int getRadius() const;
      bool requiresParameterFile() const override { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
      Variable::Type mVariable;
      int mRadius;
      Util::StatType mStatType;
      float mQuantile;
};
#endif
