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
      int    getLatDim(int iFile) const;
      int    getLonDim(int iFile) const;
      int    getTimeDim(int iFile) const;
      int    getCoefficientDim(int iFile) const;
      std::vector<int> getDims(int iFile, int iVar) const;

      int    getDim(int iFile, std::string iDim) const;
      int    getVar(int iFile, std::string iVar) const;
      int    getDimSize(int iFile, int iDim) const;
      static bool   hasDim(int iFile, std::string iDim);
      static bool   hasVar(int iFile, std::string iVar);
      std::string mDimName;
      std::string mVarName;
      void handleNetcdfError(int status, std::string message="") const;
      //! Convert linear index 'i' to vector 'iInidices'. 'iCount' specifies the size of the data
      //! Using row-major ordering (last index varies fastest)
      std::vector<int> getIndices(int i, const std::vector<int>& iCount) const;

      // Read variable from file, convert missing values
      // User must release memory
      float* getNcFloats(int iFile, int iVar);

      // Reads lat/lon/elev values from the variable, regardless of the order of the x, y dimensions
      // For lat/lon, the variable can have one or two dimensions. If the former, then the values
      // are assumed to be the same across the other dimension.
      vec2 getGridValues(int iFile, int iVariable) const;
};
#endif
