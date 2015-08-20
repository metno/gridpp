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
      enum OperatorType {
         OperatorMean      = 0,
         OperatorStd       = 30,
         OperatorQuantile  = 40
      };
      // Compute the statistic over the neighbourhood. Removes missing values.
      static float compute(const std::vector<float>& neighbourhood, OperatorType iOperator, float iQuantile=Util::MV);
   private:
      bool calibrateCore(File& iFile) const;
      Variable::Type mVariable;
      int mRadius;
      OperatorType mOperator;
      float mQuantile;
};
#endif
