#ifndef CALIBRATOR_QQ_H
#define CALIBRATOR_QQ_H
#include "Calibrator.h"
class ParameterFile;

//! Applies polynomial regression to forecasts
class CalibratorQq : public Calibrator {
   public:
      CalibratorQq(const ParameterFile* iParameterFile, Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "qq";};
   private:
      bool calibrateCore(File& iFile) const;
      const ParameterFile* mParameterFile;
      Variable::Type mVariable;
      float mLowerQuantile;
      float mUpperQuantile;
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
