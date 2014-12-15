#ifndef DATA_H
#define DATA_H
#include <netcdf.hh>
#include <vector>
#include <map>
#include "Site.h"
typedef std::vector<std::vector<std::vector<float> > > Field;

//! Represents a Netcdf data file
class DataFile {
   public:
      DataFile(std::string iFilename);
      ~DataFile();
      const Field& getField(std::string iVariable, int iTime) const;
      void setField(Field& iField, std::string iVariable, int iTime) const;
      int getNumX() const;
      int getNumY() const;
      int getNumE() const;

      //! Get the number of offsets available in the file
      int getNumOffsets() const {return mNumOffsets;};
      //! Get the initialization date (YYYYMMDD)
      int getDate() const;
      //! Get the initialization time (HH)
      int getInit() const;
      void read();
      void write();
   private:
      std::string mFilename;
      int mNumOffsets;
      NcFile mFile;
      std::vector<float> mLats;
      std::vector<float> mLons;
      std::vector<float> mElevs;
      mutable std::map<std::string, std::vector<Field> > mValues; // Variable, offset

      int mDate;
      int mInit;
};
#endif
