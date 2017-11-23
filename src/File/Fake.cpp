#include "Fake.h"
FileFake::FileFake(const Options& iOptions) :
      File("", iOptions) {
   mNLat = 10;
   mNLon = 10;
   mNEns = 2;
   mNTime = 10;

   iOptions.getValue("nLat", mNLat);
   iOptions.getValue("nLon", mNLon);
   iOptions.getValue("nEns", mNEns);
   iOptions.getValue("nTime", mNTime);
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
   mLandFractions.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      mLandFractions[i].resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         mLandFractions[i][j] = (float) i / getNumLat() / 2 + (float) j / getNumLon();
      }
   }
}

FieldPtr FileFake::getFieldCore(const Variable& iVariable, int iTime) const {
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
void FileFake::writeCore(std::vector<Variable> iVariables) {
   Util::warning("Cannot write file using the 'Fake' format");
}

bool FileFake::hasVariableCore(const Variable& iVariable) const {
   return true;
}
std::string FileFake::description() {
   std::stringstream ss;
   ss << "   type=fake                    Fake file" << std::endl;
   return ss.str();
}
