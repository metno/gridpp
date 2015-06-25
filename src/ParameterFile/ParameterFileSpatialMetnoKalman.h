#ifndef PARAMETER_FILE_SPATIAL_METNO_KALMAN_H
#define PARAMETER_FILE_SPATIAL_METNO_KALMAN_H
#include <iostream>
#include <map>
#include "ParameterFileSpatial.h"
#include "../Parameters.h"
#include "../Location.h"
class ParameterFileMetnoKalman : public ParameterFileSpatial {
   public:
      ParameterFileMetnoKalman(std::string iFilename);

      //! Get parameters for a specific time and location index
      //! @param iLocation Must be >= 0 and < getNumLocations()
      //! @return Parameters valid for time and location. If no parameters are available
      //! returns an empty parameter set
      Parameters getParameters(int iTime, const Location& iLocation) const;

      // Get unique locations in parameter set
      std::vector<Location> getLocations() const;

      // Get unique times in parameter set
      std::vector<int> getTimes() const;

      //! Set the parameter valid for specified time and location
      void setParameters(float iValue, int iTime, const Location& iLocation);

      //! Write parameter file to disk
      //! @param iFilename Write to this filename. If empty, write to the file that was read from.
      // void write(const std::string& iFilename="") const;

      //! Returns the filename where parameters are retrieved from
      std::string getFilename() const;
   private:
      std::map<Location, std::vector<float> > mParameters; // Location, Parameters (for all offsets)
      std::string mFilename;
      std::vector<int> mTimes;
      float mLocalMV;
};
#endif
