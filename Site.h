#ifndef SITE_H
#define SITE_H

//! Represents a geographic location
class Site {
   public:
      //! iId:   Station identifier
      //! iLat:  Latitude in degrees (+ for North)
      //! iLon:  Longitude in degrees (+ for East)
      //! iElev: True (i.e. not model) elevation
      Site(int iId, float iLat, float iLon, float iElev);
      float getId() const {return mId;};
      float getLat() const {return mLat;};
      float getLon() const {return mLon;};
      float getElev() const {return mElev;};
   private:
      int mId;
      float mLat;
      float mLon;
      float mElev;
};
#endif

