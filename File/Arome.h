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
      FileArome(std::string iFilename);
      ~FileArome();

      // Dimension sizes
      int getNumLat() const {return mNLat;};
      int getNumLon() const {return mNLon;};
      int getNumEns() const {return 1;};
      int getNumTime() const {return mNTime;};

      std::string getVariableName(Variable::Type iVariable) const;
      vec2 getLats() const;
      vec2 getLons() const;
      vec2 getElevs() const;
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      vec2 mLats;
      vec2 mLons;
      vec2 mElevs;

      std::vector<int> mTimes;
      vec2 getLatLonVariable(std::string iVariable) const;

      int mNTime;
      int mNLat;
      int mNLon;
};
#endif
