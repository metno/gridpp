#ifndef FILE_H
#define FILE_H
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "../Variable.h"
#include "../Uuid.h"
#include "../Field.h"

class Options;

// 3D array of data: [lat][lon][ensemble_member]
typedef std::vector<std::vector<float> > vec2; // Lat, Lon

//! Represents a data file containing spatial and temporal data initialized at one particular time
class File {
   public:
      File(std::string iFilename);
      virtual ~File();

      //! Insantiates a file. Returns null if file does not exist or cannot be parsed.
      static File* getScheme(std::string iFilename, const Options& iOptions, bool iReadOnly=false);

      FieldPtr getField(Variable::Type iVariable, int iTime) const;

      //! Get a new field initialized with missing values
      FieldPtr getEmptyField(float iFillValue=Util::MV) const;

      // Add a field to the file, overwriting existing ones (if necessary)
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
      vec2 getLandFractions() const;
      bool setLats(vec2 iLats);
      bool setLons(vec2 iLons);
      bool setElevs(vec2 iElevs);
      bool setLandFractions(vec2 iLandFractions);

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

      //! Returns a tag that uniquely identifies the latitude/longitude grid
      //! If the grid changes, a new tag is issued. Two files with the same grid
      //! will not have the same unique tag.
      Uuid getUniqueTag() const;

      //! Set the time that the file is issued
      //! @ param iTime The number of seconds since 1970-01-01 00:00:00 +00:00
      void setReferenceTime(double iTime);
      double getReferenceTime() const;
      //! Set the times of the timesteps in the file
      //! @ param iTimes vector of number of seconds since 1970-01-01 00:00:00 +00:00
      void setTimes(std::vector<double> iTimes);
      std::vector<double> getTimes() const;
      static std::string getDescriptions();
   protected:
      virtual FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const = 0;
      virtual void writeCore(std::vector<Variable::Type> iVariables) = 0;
      //! Can the subclass provide this variable?
      virtual bool hasVariableCore(Variable::Type iVariable) const = 0;

      // Subclasses must fill these fields in the constructor:
      vec2 mLats;
      vec2 mLons;
      vec2 mElevs;
      vec2 mLandFractions;
      int mNTime;
      int mNLat;
      int mNLon;
      int mNEns;
   private:
      std::string mFilename;
      mutable std::map<Variable::Type, std::vector<FieldPtr> > mFields;  // Variable, offset
      mutable Uuid mTag;
      void createNewTag() const;
      FieldPtr getEmptyField(int nLat, int nLon, int nEns, float iFillValue=Util::MV) const;
      double mReferenceTime;
      std::vector<double> mTimes;
	  static Uuid mNextTag;
};
#include "Netcdf.h"
#include "Fake.h"
#include "Point.h"
#include "NorcomQnh.h"
#endif
