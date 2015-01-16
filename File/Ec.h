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
      FileEc(std::string iFilename);
      ~FileEc();

      // Dimension sizes
      int getNumLat() const {return mNLat;};
      int getNumLon() const {return mNLon;};
      int getNumEns() const {return mNEns;};
      int getNumTime() const {return mNTime;};

      std::string getVariableName(Variable::Type iVariable) const;
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      std::vector<float> mLats;
      std::vector<float> mLons;
      std::vector<float> mElevs;

      std::vector<int> mTimes;

      int mNTime;
      int mNEns;
      int mNLat;
      int mNLon;
};
#endif
