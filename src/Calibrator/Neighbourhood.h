#ifndef CALIBRATOR_NEIGHBOURHOOD_H
#define CALIBRATOR_NEIGHBOURHOOD_H
#include "Calibrator.h"
#include "../Variable.h"
#include "../Util.h"

class ParameterFile;
class Parameters;

//! Applies an operator to a neighbourhood
class CalibratorNeighbourhood : public Calibrator {
   public:
      CalibratorNeighbourhood(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "neighbourhood";};
      int  getRadius() const;
   private:
      bool calibrateCore(File& iFile) const;
      Variable::Type mVariable;
      int mRadius;
      Util::OperatorType mOperator;
      float mQuantile;
};
#endif
