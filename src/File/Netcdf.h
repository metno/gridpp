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
      FileNetcdf(std::string iFilename, bool iReadOnly=false);
      ~FileNetcdf();

      // Does this file contain the variable?
      bool hasVariableCore(Variable::Type iVariable) const;
      virtual std::string getVariableName(Variable::Type iVariable) const = 0;
   protected:
      float getScale(NcVar* iVar) const;
      float getOffset(NcVar* iVar) const;
      NcFile mFile;

      NcDim* getDim(std::string iDim) const;
      NcVar* getVar(std::string iVar) const;
      bool   hasDim(std::string iDim) const;
      bool   hasVar(std::string iVar) const;
      static bool   hasDim(const NcFile& iFile, std::string iDim);
      static bool   hasVar(const NcFile& iFile, std::string iVar);
      static float getMissingValue(const NcVar* iVar);
};
#include "Ec.h"
#include "Arome.h"
#endif
