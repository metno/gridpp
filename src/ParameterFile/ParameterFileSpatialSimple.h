#ifndef PARAMETER_FILE_SPATIAL_SIMPLE_H
#define PARAMETER_FILE_SPATIAL_SIMPLE_H
#include <iostream>
#include <map>
#include "ParameterFileSpatial.h"
#include "../Parameters.h"
#include "../Location.h"

//! Parameters are read from a specified text file with the following format:
//! <time> <lat> <lon> <elevation> <param1> <param2> ... <paramN>
//! 0 60 10 200 3.4 2.1 5.2
//! 1 60 10 200 1.4 4.1 2.2
//! 1 58 8 150 9.2 12.2 -2.1
//! A location does not need to have parameters for all times
class ParameterFileSpatialSimple : public ParameterFileSpatial {
   public:
      ParameterFileSpatialSimple(std::string iFilename);

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
      void setParameters(Parameters iParameters, int iTime, const Location& iLocation);

      //! Write parameter file to disk
      //! @param iFilename Write to this filename. If empty, write to the file that was read from.
      void write(const std::string& iFilename="") const;

      //! Returns the filename where parameters are retrieved from
      std::string getFilename() const;
   private:
      std::map<int, std::map<Location, Parameters> > mParameters; // Offset, Location, Parameters
      std::string mFilename;
};
#endif
