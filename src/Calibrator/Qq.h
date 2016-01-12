#ifndef CALIBRATOR_QQ_H
#define CALIBRATOR_QQ_H
#include "Calibrator.h"
class ParameterFile;

//! Applies polynomial regression to forecasts
class CalibratorQq : public Calibrator {
   public:
      CalibratorQq(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "qq";};
      Parameters train(const std::vector<ObsEns>& iData) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      Variable::Type mVariable;
      float mLowerQuantile;
      float mUpperQuantile;
      std::vector<float> mQuantiles;
      struct ExtrapolationPolicy {
         enum Policy {
            OneToOne = 0,
            MeanSlope = 10,
            NearestSlope = 20,
            Zero = 30,
         };
      };
      ExtrapolationPolicy::Policy mPolicy;
      // Separate the vector of obs,fcst,obs,fcst,... into separate vectors
      static void separate(const Parameters& iParameters, std::vector<float>& iObs, std::vector<float>& iFcst);
};
#endif
