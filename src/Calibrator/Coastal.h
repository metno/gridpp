#ifndef CALIBRATOR_COASTAL_H
#define CALIBRATOR_COASTAL_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;
class Grid;

class CalibratorCoastal : public Calibrator {
   public:
      CalibratorCoastal(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "coastal";};
      Parameters train(const std::vector<ObsEnsField>& iData, const Grid& iObsGrid, const Grid& iEnsGrid, int iIobs, int iJobs, int iIEns, int iJEns) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      Variable::Type mVariable;
      int mSearchRadius;
      float mMinLsmDiff;
};
#endif
