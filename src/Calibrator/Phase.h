#ifndef CALIBRATOR_PHASE_H
#define CALIBRATOR_PHASE_H
#include "Calibrator.h"
class ParameterFile;

//! Creates a precipitation-phase field
class CalibratorPhase : public Calibrator {
   public:
      CalibratorPhase(const Options& iOptions);
      static std::string description();
      std::string name() const {return "phase";};
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
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const override;
      float mMinPrecip;
      bool mUseWetbulb;
      //! If true compute pressure based on standard atmosphere (instead of using forecasted data)
      //! This is likely a good enough approximation when computing wetbulb temperature and saves memory.
      bool mEstimatePressure;
};
#endif
