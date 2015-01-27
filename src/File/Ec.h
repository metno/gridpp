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

      // Dimension sizes
      int getNumLat() const {return mNLat;};
      int getNumLon() const {return mNLon;};
      int getNumEns() const {return mNEns;};
      int getNumTime() const {return mNTime;};

      vec2 getLats() const;
      vec2 getLons() const;
      vec2 getElevs() const;

      std::string getVariableName(Variable::Type iVariable) const;
      static bool isValid(std::string iFilename);
      std::string name() const {return "ec";};
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      vec2 mLats;
      vec2 mLons;
      vec2 mElevs;

      std::vector<int> mTimes;
      vec2 getLatLonVariable(std::string iVariable) const;

      int mNTime;
      int mNEns;
      int mNLat;
      int mNLon;
};
#endif
