#ifndef DATA_H
#define DATA_H
#include <netcdf.hh>
#include <vector>
#include <map>
#include "Variable.h"

// 3D array of data: [lat][lon][ensemble_member]
typedef std::vector<std::vector<std::vector<float> > > Field;

//! Represents a Netcdf data file
class DataFile {
   public:
      DataFile(std::string iFilename);
      ~DataFile();
      const Field& getField(Variable::Type iVariable, int iTime) const;

      // Schedule field to be written
      void   addField(Field& iField, Variable::Type iVariable, int iTime);
      void write();

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
      std::string getVariableName(Variable::Type iVariable) const;
      // Does this file contain the variable?
      bool hasVariable(Variable::Type iVariable) const;

   private:
      void loadFields(Variable::Type iVariable) const;
      void saveField(Field* iField, Variable::Type iVariable, int iTime) const;
      mutable std::map<Variable::Type, std::vector<Field*> > mReadFields;  // Variable, offset

      mutable std::map<Variable::Type, std::vector<Field*> > mWriteFields; // Variable, offset
      float getScale(NcVar* iVar) const;
      float getOffset(NcVar* iVar) const;
      std::string mFilename;
      int mNumOffsets;
      NcFile mFile;
      std::vector<float> mLats;
      std::vector<float> mLons;
      std::vector<float> mElevs;
      // Creates empty field of 0s

      int mNTime;
      int mNEns;
      int mNLat;
      int mNLon;

      int mDate;
      int mInit;
};
#endif
