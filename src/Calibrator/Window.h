#ifndef CALIBRATOR_WINDOW_H
#define CALIBRATOR_WINDOW_H
#include "Calibrator.h"
#include "../Variable.h"

// Applies a statistical operator to values within a temporal window
class CalibratorWindow : public Calibrator {
   public:
      CalibratorWindow(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const override {return "window";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
      Variable::Type mVariable;
      int mRadius;
      Util::StatType mStatType;
      float mQuantile;
};
#endif
