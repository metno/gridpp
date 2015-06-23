#ifndef LOCATION_H
#define LOCATION_H
#include "Util.h"

//! Represents a point in space
class Location {
   public:
      Location(float iLat, float iLon, float iElev=Util::MV);
      //! Accessors
      float lat() const;
      float lon() const;
      float elev() const;
      //! Mutators
      void lat(float iLat);
      void lon(float iLon);
      void elev(float iElev);
      //! Used for sorting in std::map
      bool operator<(const Location &right) const;
      bool operator==(const Location &right) const;
   private:
      float mLat;
      float mLon;
      float mElev;
};
#endif
