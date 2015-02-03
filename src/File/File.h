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
      int getNumLat() const;
      int getNumLon() const;
      int getNumEns() const;
      int getNumTime() const;
      vec2 getLats() const;
      vec2 getLons() const;
      vec2 getElevs() const;
      bool setLats(vec2 iLats);
      bool setLons(vec2 iLons);
      bool setElevs(vec2 iElevs);

      //! Does this file provide the variable (deriving it if necessary)?
      bool hasVariable(Variable::Type iVariable) const;

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
      virtual FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const = 0;
      virtual void writeCore(std::vector<Variable::Type> iVariables) = 0;
      //! Can the subclass provide this variable?
      virtual bool hasVariableCore(Variable::Type iVariable) const = 0;

      // Subclasses must fill these fields in the constructor:
      vec2 mLats;
      vec2 mLons;
      vec2 mElevs;
      int mNTime;
      int mNLat;
      int mNLon;
      int mNEns;
   private:
      std::string mFilename;
      mutable std::map<Variable::Type, std::vector<FieldPtr> > mFields;  // Variable, offset
      mutable boost::uuids::uuid mTag;
      void createNewTag() const;
      FieldPtr getEmptyField(int nLat, int nLon, int nEns, float iFillValue=Util::MV) const;
};
#include "Netcdf.h"
#include "Fake.h"
#endif
