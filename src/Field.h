#ifndef FIELD_H
#define FIELD_H
#include <boost/shared_ptr.hpp>
#include <vector>
#include "Util.h"

//! Encapsulates gridded data in 3 dimensions: latitude, longitude, ensemble member
class Field {
   public:
      //! Initialize 3D field
      //! @param nLat number of latitudes
      //! @param nLon number of longitudes
      //! @param nLon number of ensemble members
      //! @param iFillValue initialize all values in field with this
      Field(int nLat, int nLon, int nEns, float iFillValue=Util::MV);
      //! Access to data
      //! @param i latitude index
      //! @param j longitude index
      //! @param k ensemble index
      //! @return data at specified coordinate
      float      & operator()(unsigned int i, unsigned int j, unsigned int k);
      float const& operator()(unsigned int i, unsigned int j, unsigned int k) const;
      std::vector<float> operator()(unsigned int i, unsigned int j) const;
      //! Are all values (for all lat/lon/ens) in fields identical?
      bool operator==(const Field& iField) const;
      bool operator!=(const Field& iField) const;
      int getNumLat() const;
      int getNumLon() const;
      int getNumEns() const;
   protected:
      //! Data values. Index for ensemble changes fastest.
      std::vector<float> mValues;
      int mNLat;
      int mNLon;
      int mNEns;
      //! Index into flat array that corresponds to coordinate
      int getIndex(unsigned int i, unsigned int j, unsigned int k) const;
};
typedef boost::shared_ptr<Field> FieldPtr;
#endif
