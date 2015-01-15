#ifndef FILE_H
#define FILE_H
#include <netcdf.hh>
#include <vector>
#include <map>
#include "Variable.h"

// 3D array of data: [lat][lon][ensemble_member]
typedef std::vector<std::vector<std::vector<float> > > Field;

//! Represents a Netcdf data file
class File {
   public:
      File(std::string iFilename);
      ~File();
      Field& getField(Variable::Type iVariable, int iTime) const;

      // Add an auxillary field that user of the class can read.
      void addField(Field& iField, Variable::Type iVariable, int iTime);

      // Write these variables to file
      void write(std::vector<Variable::Type> iVariables);

      Field& getEmptyField(int nLat, int nLon, int nEns) const;
      Field& getEmptyField() const;

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
      void loadFields(Variable::Type iVariable) const;
      void saveField(Field* iField, Variable::Type iVariable, int iTime) const;
      mutable std::map<Variable::Type, std::vector<Field*> > mFields;  // Variable, offset
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
