#ifndef FILE_EC_H
#define FILE_EC_H
#include <netcdf.hh>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "Netcdf.h"
#include "../Variable.h"

//! Represents a Netcdf data file
class FileEc : public FileNetcdf {
   public:
      FileEc(std::string iFilename, bool iReadOnly=false);

      std::string getVariableName(Variable::Type iVariable) const;
      static bool isValid(std::string iFilename);
      std::string name() const {return "ec";};
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;

      std::vector<int> mTimes;
      vec2 getLatLonVariable(std::string iVariable) const;
};
#endif
