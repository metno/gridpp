#include "Fake.h"
FileFake::FileFake(int nLat, int nLon, int nEns, int nTime) :
      File("") {
   mNLat = nLat;
   mNLon = nLon;
   mNEns = nEns;
   mNTime = nTime;
   if(!Util::isValid(mNLat) || mNLat <= 0) {
      Util::error("FileFake: Invalid number of latitudes");
   }
   if(!Util::isValid(mNLon) || mNLon <= 0) {
      Util::error("FileFake: Invalid number of longitudes");
   }
   if(!Util::isValid(mNEns) || mNEns <= 0) {
      Util::error("FileFake: Invalid number of ensemble members");
   }
   if(!Util::isValid(mNTime) || mNTime <= 0) {
      Util::error("FileFake: Invalid number of times");
   }
   mLats.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      mLats[i].resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         mLats[i][j] = 50 + 10.0 * i / getNumLat();
      }
   }
   mLons.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      mLons[i].resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         mLons[i][j] = 0 + 10.0 * j / getNumLon();
      }
   }
   mElevs.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      mElevs[i].resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         mElevs[i][j] = 0;
      }
   }
}

FieldPtr FileFake::getFieldCore(Variable::Type iVariable, int iTime) const {

   FieldPtr field = getEmptyField();

   for(int i = 0; i < getNumLat(); i++) {
      for(int j = 0; j < getNumLon(); j++) {
         for(int e = 0; e < getNumEns(); e++) {
            (*field)(i,j,e) = i + j + e;
         }
      }
   }
   return field;
}
void FileFake::writeCore(std::vector<Variable::Type> iVariables) {
   Util::warning("Cannot write file using the 'Fake' format");
}

bool FileFake::hasVariableCore(Variable::Type iVariable) const {
   if(iVariable == Variable::PrecipAcc)
      return false;
   else
      return true;
}
