#include "Location.h"

Location::Location(float iLat, float iLon, float iElev) :
      mLat(iLat),
      mLon(iLon),
      mElev(iElev) {

}

bool Location::operator<(const Location &right) const {
   if(mLat == right.lat()) {
      if(mLon == right.lon()) {
         return mElev < right.mElev;
      }
      else {
         return mLon < right.mLon;
      }
   }
   else {
      return mLat < right.mLat;
   }
}

bool Location::operator==(const Location &right) const {
   return (mLat == right.lat() && mLon == right.lon() && mElev == right.elev());
}
float Location::lat() const {
   return mLat;
}
float Location::lon() const {
   return mLon;
}
float Location::elev() const {
   return mElev;
}
void Location::lat(float iLat) {
   mLat = iLat;
}
void Location::lon(float iLon) {
   mLon = iLon;
}
void Location::elev(float iElev) {
   mElev = iElev;
}
