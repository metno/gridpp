#ifndef FILE_EC_H
#define FILE_EC_H
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "NetcdfBase.h"
#include "../Variable.h"
#include "../Options.h"

//! Represents an ensemble Netcdf data file from ECMWF
//! Must have:
//!    Time dim: time
//!    Latitude dim:  lat, latitude, y
//!    Longitude dim: lon, longitude, x
//!    Latitude var:  lat, latitude, y. With the latitude dim, and possibly the lon dim
//!    Longitude var: lon, longitude, x. With the longitude dim, and possibly the lat dim
//!    5 dimensional variables (time, *, ensemble_member, <lat dim>, <lon dim>)
class FileEc : public FileNetcdfBase {
   public:
      FileEc(std::string iFilename, const Options& iOptions=Options(), bool iReadOnly=false);

      std::string getVariableName(Variable::Type iVariable) const;
      static bool isValid(std::string iFilename);
      static std::string description();
      std::string name() const {return "ec";};
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      FieldPtr getFieldCore(std::string iVariable, int iTime) const;

      std::vector<int> mTimes;
      vec2 getGridValues(int iVariable) const;
      void writeAltitude() const;
      void defineAltitude();
      int getLatDim() const;
      int getLonDim() const;
      int getLatVar() const;
      int getLonVar() const;
};
#endif
