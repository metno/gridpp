#ifndef CALIBRATOR_PHASE_H
#define CALIBRATOR_PHASE_H
#include "Calibrator.h"
class ParameterFile;

// Ensure that if a member has precip it also has full cloud cover
class CalibratorPhase : public Calibrator {
   public:
      CalibratorPhase(const ParameterFile* iParameterFile);
      static std::string description();
      std::string name() const {return "phase";};
      enum Phase {PhaseNone = 0, PhaseRain = 1, PhaseSleet = 2, PhaseSnow = 3};
      //! Compute wetbulb temperature
      //! @param iTemperature Temperature in K
      //! @param iPressure Pressure in pa
      //! @param iRelativeHumidity Relative humidity (out of 1)
      //! @return Wetbulb temperature in K
      static float getWetbulb(float iTemperature, float iPressure, float iRelativeHumidity);

      float getMinPrecip() const;
      void  setMinPrecip(float iMinPrecip);
      void  setUseWetbulb(bool iUseWetbulb);
      bool  getUseWetbulb();
   private:
      bool calibrateCore(File& iFile) const;
      const ParameterFile* mParameterFile;
      float mMinPrecip;
      bool mUseWetbulb;
};
#endif
