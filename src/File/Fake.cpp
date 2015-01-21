#include "Fake.h"
FileFake::FileFake(int nLat, int nLon, int nEns, int nTime) :
      File(""), mNLat(nLat), mNLon(nLon), mNEns(nEns), mNTime(nTime) {
}

FieldPtr FileFake::getFieldCore(Variable::Type iVariable, int iTime) const {

   FieldPtr field = getEmptyField();

   for(int i = 0; i < getNumLat(); i++) {
      for(int j = 0; j < getNumLon(); j++) {
         for(int e = 0; e < getNumEns(); e++) {
            (*field)[i][j][e] = i + j + e;
         }
      }
   }
   return field;
}

vec2 FileFake::getLats() const {
   vec2 lats;
   lats.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      lats.resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         lats[i][j] = 50 + 10.0 * i / getNumLat();
      }
   }
   return lats;
}
vec2 FileFake::getLons() const {
   vec2 lons;
   lons.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      lons.resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         lons[i][j] = 0 + 10.0 * j / getNumLon();
      }
   }
   return lons;
}
vec2 FileFake::getElevs() const {
   vec2 elevs;
   elevs.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      elevs.resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         elevs[i][j] = 0;
      }
   }
   return elevs;
}

bool FileFake::setLats(vec2 iLats) {
   if(iLats.size() != mNLat || iLats[0].size() != mNLon)
      return false;
   mLats = iLats;
   return true;
}
bool FileFake::setLons(vec2 iLons) {
   if(iLons.size() != mNLat || iLons[0].size() != mNLon)
      return false;
   mLons = iLons;
   return true;
}
bool FileFake::setElevs(vec2 iElevs) {
   if(iElevs.size() != mNLat || iElevs[0].size() != mNLon)
      return false;
   mElevs = iElevs;
   return true;
}
