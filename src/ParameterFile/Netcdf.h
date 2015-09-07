#ifndef PARAMETER_FILE_NETCDF_H
#define PARAMETER_FILE_NETCDF_H
#include <iostream>
#include <netcdf.h>
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

      void write() const;
   private:
      std::vector<int> mTimes;
      float mLocalMV;
      int    getDim(int iFile, std::string iDim) const;
      int    getVar(int iFile, std::string iVar) const;
      int    getDimSize(int iFile, int iDim) const;
      static bool   hasDim(int iFile, std::string iDim);
      static bool   hasVar(int iFile, std::string iVar);
      std::string mDimName;
      std::string mVarName;
      void handleNetcdfError(int status, std::string message="") const;
};
#endif
