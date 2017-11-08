#ifndef FILE_AROME_H
#define FILE_AROME_H
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "NetcdfBase.h"
#include "../Variable.h"
#include "../Options.h"

//! Represents a Netcdf data file
class FileArome : public FileNetcdfBase {
   public:
      FileArome(std::string iFilename, const Options& iOptions=Options(), bool iReadOnly=false);
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
      // Must be one of "latitude", "longitude", or "altitude"
      vec2 getLatLonVariable(std::string iVariable) const;
      void writeLatLonVariable(std::string iVariable);
      void defineLatLonVariable(std::string iVariable);
      int mDate;
      // What model level does this variable come from?
      int getLevel(std::string iVariable) const;
      std::string mXName;
      std::string mYName;
      std::string mLatName;
      std::string mLonName;
};
#endif
