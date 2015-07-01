#ifndef PARAMETER_FILE_NETCDF_H
#define PARAMETER_FILE_NETCDF_H
#include <iostream>
#include <netcdf.hh>
#include <map>
#include "ParameterFile.h"
#include "../Parameters.h"
#include "../Location.h"
#include "../File/File.h"
class ParameterFileNetcdf : public ParameterFile {
   public:
      ParameterFileNetcdf(const Options& iOptions);

      bool isFixedSize() const {return true;};
      std::vector<int> getTimes() const;

      static bool isValid(std::string iFilename);

      static std::string description();
      std::string name() const {return "netcdf";};
   private:
      std::vector<int> mTimes;
      float mLocalMV;
      NcDim* getDim(const NcFile& iFile, std::string iDim) const;
      NcVar* getVar(const NcFile& iFile, std::string iVar) const;
      bool hasDim(const NcFile& iFile, std::string iDim) const;
      bool hasVar(const NcFile& iFile, std::string iVar) const;
      std::string mDimName;
      std::string mVarName;
};
#endif
