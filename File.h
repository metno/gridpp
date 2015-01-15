#ifndef FILE_H
#define FILE_H
#include <netcdf.hh>
#include <vector>
#include <map>
#include "Variable.h"
#include <boost/shared_ptr.hpp>

// 3D array of data: [lat][lon][ensemble_member]
typedef std::vector<std::vector<std::vector<float> > > Field;
typedef boost::shared_ptr<Field> FieldPtr;

//! Represents a Netcdf data file
class File {
   public:
      File(std::string iFilename);
      ~File();

      // Accessors
      FieldPtr getField(Variable::Type iVariable, int iTime) const;
      //! Get a new field initialized with missing values
      FieldPtr getEmptyField() const;

      // Add a field to the file, overwriting existing ones (if available)
      void addField(FieldPtr iField, Variable::Type iVariable, int iTime);

      // Write these variables to file
      void write(std::vector<Variable::Type> iVariables);

      // Dimension sizes
      int getNumLat() const {return mNLat;};
      int getNumLon() const {return mNLon;};
      int getNumEns() const {return mNEns;};
      int getNumTime() const {return mNTime;};

      //! Get the initialization date (YYYYMMDD)
      int getDate() const;
      //! Get the initialization time (HH)
      int getInit() const;
      std::string getFilename() const;
      std::string getVariableName(Variable::Type iVariable) const;

      // Does this file contain the variable?
      bool hasVariable(Variable::Type iVariable) const;
      bool hasSameDimensions(const File& iOther) const;
      std::string getDimenionString() const;
   private:
      FieldPtr getEmptyField(int nLat, int nLon, int nEns) const;
      void loadFields(Variable::Type iVariable) const;
      void saveField(FieldPtr iField, Variable::Type iVariable, int iTime) const;
      mutable std::map<Variable::Type, std::vector<FieldPtr> > mFields;  // Variable, offset
      float getScale(NcVar* iVar) const;
      float getOffset(NcVar* iVar) const;
      std::string mFilename;
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

      int mDate;
      int mInit;
};
#endif
