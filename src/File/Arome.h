#ifndef FILE_AROME_H
#define FILE_AROME_H
#include <netcdf.hh>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "Netcdf.h"
#include "../Variable.h"

//! Represents a Netcdf data file
class FileArome : public FileNetcdf {
   public:
      FileArome(std::string iFilename, bool iReadOnly=false);
      ~FileArome();

      std::string getVariableName(Variable::Type iVariable) const;
      int  getDate() const;
      // Is the file readable in this format?
      static bool isValid(std::string iFilename);
      static std::string description();
      std::string name() const {return "arome";};
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      FieldPtr getFieldCore(std::string iVariable, int iTime) const;
      vec2 getLatLonVariable(std::string iVariable) const;
      int mDate;
};
#endif
