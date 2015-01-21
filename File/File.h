#ifndef FILE_H
#define FILE_H
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "../Variable.h"

// 3D array of data: [lat][lon][ensemble_member]
typedef std::vector<std::vector<std::vector<float> > > Field;
typedef std::vector<std::vector<float> > vec2; // Lat, Lon
typedef boost::shared_ptr<Field> FieldPtr;

//! Represents a Netcdf data file
class File {
   public:
      File(std::string iFilename);
      ~File();
      static File* getScheme(std::string iFilename);

      FieldPtr getField(Variable::Type iVariable, int iTime) const;
      //! Get a new field initialized with missing values
      FieldPtr getEmptyField() const;

      // Add a field to the file, overwriting existing ones (if available)
      void addField(FieldPtr iField, Variable::Type iVariable, int iTime) const;

      // Write these variables to file
      void write(std::vector<Variable::Type> iVariables);

      // Dimension sizes
      virtual int getNumLat() const  = 0;
      virtual int getNumLon() const  = 0;
      virtual int getNumEns() const  = 0;
      virtual int getNumTime() const = 0;
      virtual bool hasVariable(Variable::Type iVariable) const = 0;

      virtual vec2 getLats() const = 0;
      virtual vec2 getLons() const = 0;
      virtual vec2 getElevs() const = 0;

      std::string getFilename() const;

      bool hasSameDimensions(const File& iOther) const;
      std::string getDimenionString() const;
      void initNewVariable(Variable::Type iVariable);
   protected:
      std::string mFilename;
      virtual FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const = 0;
      virtual void writeCore(std::vector<Variable::Type> iVariables) = 0;
      FieldPtr getEmptyField(int nLat, int nLon, int nEns) const;
   private:
      mutable std::map<Variable::Type, std::vector<FieldPtr> > mFields;  // Variable, offset
};
#include "Netcdf.h"
#include "Fake.h"
#endif
