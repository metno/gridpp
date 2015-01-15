#ifndef FILE_NETCDF_H
#define FILE_NETCDF_H
#include <netcdf.hh>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "File.h"
#include "../Variable.h"

// 3D array of data: [lat][lon][ensemble_member]
typedef std::vector<std::vector<std::vector<float> > > Field;
typedef boost::shared_ptr<Field> FieldPtr;

//! Represents a Netcdf data file
class FileNetcdf : public File {
   public:
      FileNetcdf(std::string iFilename);
      ~FileNetcdf();

      // Dimension sizes
      int getNumLat() const {return mNLat;};
      int getNumLon() const {return mNLon;};
      int getNumEns() const {return mNEns;};
      int getNumTime() const {return mNTime;};

      std::string getVariableName(Variable::Type iVariable) const;

      // Does this file contain the variable?
      bool hasVariable(Variable::Type iVariable) const;
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      void loadFields(Variable::Type iVariable) const;
      float getScale(NcVar* iVar) const;
      float getOffset(NcVar* iVar) const;
      int mNumOffsets;
      NcFile mFile;
      std::vector<float> mLats;
      std::vector<float> mLons;
      std::vector<float> mElevs;
      // Creates empty field of 0s

      NcDim* getDim(std::string iDim) const;
      NcVar* getVar(std::string iVar) const;
      static float getMissingValue(const NcVar* iVar);
      std::vector<int> mTimes;

      int mNTime;
      int mNEns;
      int mNLat;
      int mNLon;
};
#endif
