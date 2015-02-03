#ifndef FILE_H
#define FILE_H
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include "../Variable.h"
#include "../Util.h"
#include "../Field.h"

// 3D array of data: [lat][lon][ensemble_member]
typedef std::vector<std::vector<float> > vec2; // Lat, Lon

//! Represents a Netcdf data file
class File {
   public:
      File(std::string iFilename);
      virtual ~File();
      static File* getScheme(std::string iFilename, bool iReadOnly=false);

      FieldPtr getField(Variable::Type iVariable, int iTime) const;
      //! Get a new field initialized with missing values
      FieldPtr getEmptyField(float iFillValue=Util::MV) const;

      // Add a field to the file, overwriting existing ones (if available)
      void addField(FieldPtr iField, Variable::Type iVariable, int iTime) const;

      // Write these variables to file
      void write(std::vector<Variable::Type> iVariables);

      // Dimension sizes
      virtual int getNumLat() const  = 0;
      virtual int getNumLon() const  = 0;
      virtual int getNumEns() const  = 0;
      virtual int getNumTime() const = 0;

      //! Does this file provide the variable (deriving it if necessary)?
      bool hasVariable(Variable::Type iVariable) const;

      virtual vec2 getLats() const = 0;
      virtual vec2 getLons() const = 0;
      virtual vec2 getElevs() const = 0;

      std::string getFilename() const;

      bool hasSameDimensions(const File& iOther) const;
      std::string getDimenionString() const;
      void initNewVariable(Variable::Type iVariable);
      virtual std::string name() const = 0;
      //! Clear the retrieved/computed fields stored in cache
      void clear();
      //! How many bytes of retrieved/computed  data are stored in cache?
      //! @return Number of bytes
      long getCacheSize() const;

      boost::uuids::uuid getUniqueTag() const;
   protected:
      std::string mFilename;
      virtual FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const = 0;
      virtual void writeCore(std::vector<Variable::Type> iVariables) = 0;
      //! Can the subclass provide this variable?
      virtual bool hasVariableCore(Variable::Type iVariable) const = 0;
      FieldPtr getEmptyField(int nLat, int nLon, int nEns, float iFillValue=Util::MV) const;
   private:
      mutable std::map<Variable::Type, std::vector<FieldPtr> > mFields;  // Variable, offset
      mutable boost::uuids::uuid mTag;
};
#include "Netcdf.h"
#include "Fake.h"
#endif
