#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H
#include <iostream>
#include <map>
#include "../Parameters.h"
#include "../Location.h"
#include "../Options.h"
#include "../Scheme.h"

//! Represents a collection of parameters, one set for each forecast time
class ParameterFile : public Scheme {
   public:
      ParameterFile(const Options& iOptions);

      //! Get the parameter valid for specified forecast timestep. This is an index, not an hour.
      Parameters getParameters(int iTime, const Location& iLocation) const;
      //! Only use this if isLocationDependent() is false
      Parameters getParameters(int iTime) const;

      static ParameterFile* getScheme(std::string iName, const Options& iOptions=Options());

      //! Does this file provide different parameters for different locations?
      virtual bool isLocationDependent() const;
      //! Does this file provide different parameters for different lead times?
      bool isTimeDependent() const;
      //! Does this file have the same number of parameters for all locations/lead times?
      virtual bool isFixedSize() const = 0;

      //! Set the parameter valid for specified time
      void setParameters(Parameters iParameters, int iTime, const Location& iLocation);
      void setParameters(Parameters iParameters, int iTime);

      std::vector<Location> getLocations() const;
      virtual std::vector<int> getTimes() const = 0;
      virtual bool isReadable() const = 0;

      // Returns Util::MV if the number are not consistent
      int getNumParameters() const;

      std::string getFilename() const;
      virtual std::string name() const = 0;
      static std::string getDescription(bool iSpatialOnly=false);
      virtual void write() const {};

      static std::string getDescriptions();
   protected:

      // Store all location-dependent parameters here
      std::map<Location, std::map<int, Parameters>, Location::CmpIgnoreElevation > mParameters; // Location, Offset, Parameters
      std::string mFilename;
      void setFilename(std::string iFilename);
};
#include "MetnoKalman.h"
#include "Text.h"
#include "Simple.h"
#include "Netcdf.h"
#endif
