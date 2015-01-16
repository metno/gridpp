#ifndef FILE_NETCDF_H
#define FILE_NETCDF_H
#include <netcdf.hh>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "File.h"
#include "../Variable.h"

//! Represents a Netcdf data file
class FileNetcdf : public File {
   public:
      FileNetcdf(std::string iFilename);
      ~FileNetcdf();

      // Does this file contain the variable?
      bool hasVariable(Variable::Type iVariable) const;
      virtual std::string getVariableName(Variable::Type iVariable) const = 0;
   protected:
      float getScale(NcVar* iVar) const;
      float getOffset(NcVar* iVar) const;
      NcFile mFile;

      NcDim* getDim(std::string iDim) const;
      NcVar* getVar(std::string iVar) const;
      static float getMissingValue(const NcVar* iVar);
};
#include "Ec.h"
#include "Arome.h"
#endif
