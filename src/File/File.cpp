#include "File.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include "../Util.h"
#include "../Options.h"
Uuid File::mNextTag = 0;

File::File(std::string iFilename) :
      mFilename(iFilename),
      mReferenceTime(Util::MV) {
   createNewTag();
}

File* File::getScheme(std::string iFilename, const Options& iOptions, bool iReadOnly) {
   File* file;
   // Determine the filetype, either through user-specified option type=...
   // or by autodetecting.
   std::string type = "";
   if(!iOptions.getValue("type", type)) {
      // Autodetect type based on content
      if(FileArome::isValid(iFilename)) {
         type = "arome";
      }
      else if(FileEc::isValid(iFilename)) {
         type = "ec";
      }
   }

   // Instantiate the file
   if(type == "") {
      if(!Util::exists(iFilename)) {
         Util::warning("File '" + iFilename + " does not exist");
         return NULL;
      }
      else {
         Util::warning("Could not find suitable parser for '" + iFilename + "'");
         return NULL;
      }
   }
   else if(type == "arome") {
      file = new FileArome(iFilename, iReadOnly);
   }
   else if(type == "ec") {
      file = new FileEc(iFilename, iReadOnly);
   }
   else if(type == "point") {
      file = new FilePoint(iFilename, iOptions);
   }
   else if(type == "text") {
      file = new FileText(iFilename, iOptions);
   }
   else if(type == "norcomQnh") {
      file = new FileNorcomQnh(iFilename, iOptions);
   }
   else {
      Util::error("Could not understand file type " + type);
   }
   return file;
}

FieldPtr File::getField(Variable::Type iVariable, int iTime) const {
   // Determine if values have been cached
   std::map<Variable::Type, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
   bool needsReading = it == mFields.end();
   if(!needsReading) {
      if(mFields[iVariable].size() <= iTime) {
         std::stringstream ss;
         ss << "Attempted to access variable '" << Variable::getTypeName(iVariable) << "' for time " << iTime
            << " in file '" << getFilename() << "'";
         Util::error(ss.str());
      }

      needsReading = mFields[iVariable][iTime] == NULL;
   }
   else {
      mFields[iVariable].resize(getNumTime());
   }

   if(needsReading) {
      // Load non-derived variable from file
      if(hasVariableCore(iVariable)) {
         for(int t = 0; t < getNumTime(); t++) {
            mFields[iVariable][t] = getFieldCore(iVariable, t);
         }
      }
      // Try to derive the field
      else if(iVariable == Variable::Precip) {
         if(!hasVariableCore(Variable::PrecipAcc)) {
            Util::error("Cannot derive Precip");
         }
         // Deaccumulate
         FieldPtr field = getEmptyField();
         addField(field, Variable::Precip, 0); // First offset is 0

         for(int t = 1; t < getNumTime(); t++) {
            FieldPtr field = getEmptyField();
            const FieldPtr acc0  = getField(Variable::PrecipAcc, t-1);
            const FieldPtr acc1  = getField(Variable::PrecipAcc, t);
            for(int lat = 0; lat < getNumLat(); lat++) {
               for(int lon = 0; lon < getNumLon(); lon++) {
                  for(int e = 0; e < getNumEns(); e++) {
                     float a1 = (*acc1)(lat,lon,e);
                     float a0 = (*acc0)(lat,lon,e);
                     float value = Util::MV;
                     if(Util::isValid(a1) && Util::isValid(a0)) {
                         value = a1 - a0;
                         if(value < 0)
                            value = 0;
                     }
                     (*field)(lat,lon,e) = value;
                  }
               }
            }
            addField(field, Variable::Precip, t);
         }
      }
      else if(iVariable == Variable::PrecipAcc) {
         if(!hasVariableCore(Variable::Precip)) {
            Util::error("Cannot derive PrecipAcc");
         }
         // Accumulate
         FieldPtr prevAccum = getEmptyField(0);
         addField(prevAccum, Variable::PrecipAcc, 0); // First offset is 0

         for(int t = 1; t < getNumTime(); t++) {
            FieldPtr currAccum = getEmptyField();
            const FieldPtr currPrecip  = getField(Variable::Precip, t);
            for(int lat = 0; lat < getNumLat(); lat++) {
               for(int lon = 0; lon < getNumLon(); lon++) {
                  for(int e = 0; e < getNumEns(); e++) {
                     float a = (*prevAccum)(lat,lon,e);
                     float p = (*currPrecip)(lat,lon,e);
                     float value = Util::MV;
                     if(Util::isValid(a) && Util::isValid(p)) {
                         value = a + p;
                         if(value < 0)
                            value = 0;
                     }
                     (*currAccum)(lat,lon,e) = value;
                  }
               }
            }
            addField(currAccum, Variable::PrecipAcc, t);
            prevAccum = currAccum;
         }
      }
      else if(iVariable == Variable::W) {
         if(hasVariableCore(Variable::U) && hasVariableCore(Variable::V)) {
            for(int t = 0; t < getNumTime(); t++) {
               FieldPtr windSpeed = getEmptyField();
               const FieldPtr u = getField(Variable::U, t);
               const FieldPtr v = getField(Variable::V, t);
               for(int lat = 0; lat < getNumLat(); lat++) {
                  for(int lon = 0; lon < getNumLon(); lon++) {
                     for(int e = 0; e < getNumEns(); e++) {
                        float currU = (*u)(lat,lon,e);
                        float currV = (*v)(lat,lon,e);
                        (*windSpeed)(lat,lon,e) = sqrt(currU*currU + currV*currV);
                     }
                  }
               }
               addField(windSpeed, Variable::W, t);
            }
         }
         else if(hasVariableCore(Variable::Xwind) && hasVariableCore(Variable::Ywind)) {
            for(int t = 0; t < getNumTime(); t++) {
               FieldPtr windSpeed = getEmptyField();
               const FieldPtr x = getField(Variable::Xwind, t);
               const FieldPtr y = getField(Variable::Ywind, t);
               for(int lat = 0; lat < getNumLat(); lat++) {
                  for(int lon = 0; lon < getNumLon(); lon++) {
                     for(int e = 0; e < getNumEns(); e++) {
                        float currX = (*x)(lat,lon,e);
                        float currY = (*y)(lat,lon,e);
                        (*windSpeed)(lat,lon,e) = sqrt(currX*currX + currY*currY);
                     }
                  }
               }
               addField(windSpeed, Variable::W, t);
            }
         }
         else {
            Util::error("Cannot derive wind speed from variables in file");
         }
      }
      else if(iVariable == Variable::WD) {
         if(hasVariableCore(Variable::U) && hasVariableCore(Variable::V)) {
            for(int t = 0; t < getNumTime(); t++) {
               FieldPtr windDir = getEmptyField();
               const FieldPtr u = getField(Variable::U, t);
               const FieldPtr v = getField(Variable::V, t);
               for(int lat = 0; lat < getNumLat(); lat++) {
                  for(int lon = 0; lon < getNumLon(); lon++) {
                     for(int e = 0; e < getNumEns(); e++) {
                        float currU = (*u)(lat,lon,e);
                        float currV = (*v)(lat,lon,e);
                        float dir = std::atan2(-currU,-currV) * 180 / Util::pi;
                        if(dir < 0)
                           dir += 360;
                        (*windDir)(lat,lon,e) = dir;
                     }
                  }
               }
               addField(windDir, Variable::WD, t);
            }
         }
         else if(hasVariableCore(Variable::Xwind) && hasVariableCore(Variable::Ywind)) {
            for(int t = 0; t < getNumTime(); t++) {
               FieldPtr windDir = getEmptyField();
               const FieldPtr x = getField(Variable::Xwind, t);
               const FieldPtr y = getField(Variable::Ywind, t);
               for(int lat = 0; lat < getNumLat(); lat++) {
                  for(int lon = 0; lon < getNumLon(); lon++) {
                     for(int e = 0; e < getNumEns(); e++) {
                        float currX = (*x)(lat,lon,e);
                        float currY = (*y)(lat,lon,e);
                        float dir = std::atan2(-currX,-currY) * 180 / Util::pi;
                        if(dir < 0)
                           dir += 360;
                        (*windDir)(lat,lon,e) = dir;
                     }
                  }
               }
               addField(windDir, Variable::WD, t);
            }
         }
         else {
            Util::warning("Cannot derive wind speed from variables in file");
            for(int t = 0; t < getNumTime(); t++) {
               FieldPtr field = getEmptyField();
               addField(field, iVariable, t);
            }
         }
      }
      else {
         std::string variableType = Variable::getTypeName(iVariable);
         Util::warning(variableType + " not available in '" + getFilename() + "'");
         for(int t = 0; t < getNumTime(); t++) {
            FieldPtr field = getEmptyField();
            addField(field, iVariable, t);
         }
      }
   }
   if(mFields[iVariable].size() <= iTime) {
      std::stringstream ss;
      ss << "Attempted to access variable '" << Variable::getTypeName(iVariable) << "' for time " << iTime
         << " in file '" << getFilename() << "'";
      Util::error(ss.str());
   }
   FieldPtr field = mFields[iVariable][iTime];
   return field;
}

File::~File() {
}

void File::write(std::vector<Variable::Type> iVariables) {
   writeCore(iVariables);
   // mCache.clear();
}


FieldPtr File::getEmptyField(float iFillValue) const {
   return getEmptyField(getNumLat(), getNumLon(), getNumEns(), iFillValue);
}
FieldPtr File::getEmptyField(int nLat, int nLon, int nEns, float iFillValue) const {
   FieldPtr field = FieldPtr(new Field(nLat, nLon, nEns, iFillValue));
   return field;
}

void File::addField(FieldPtr iField, Variable::Type iVariable, int iTime) const {
   std::map<Variable::Type, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
   if(it == mFields.end()) {
      mFields[iVariable].resize(getNumTime());
   }

   mFields[iVariable][iTime] = iField;
}

bool File::hasSameDimensions(const File& iOther) const {
   if(getNumLat() == iOther.getNumLat()
         && getNumLon() == iOther.getNumLon()
         && getNumEns() == iOther.getNumEns()
         && getNumTime() == iOther.getNumTime())
      return true;
   return false;
}

std::string File::getDimenionString() const {
   std::stringstream ss;
   ss << "[" << getNumTime() << " " << getNumEns() << " " << getNumLat() << " " << getNumLon()<< "]";
   return ss.str();
}

std::string File::getFilename() const {
   return mFilename;
}

void File::initNewVariable(Variable::Type iVariable) {
   if(!hasVariable(iVariable)) {
      for(int t = 0; t < getNumTime(); t++) {
         addField(getEmptyField(), iVariable, t);
      }
   }
}
bool File::hasVariable(Variable::Type iVariable) const {
   bool status = hasVariableCore(iVariable);
   if(status)
      return true;
// Check if field is derivable
   if(iVariable == Variable::Precip) {
      return hasVariableCore(Variable::PrecipAcc);
   }
   else if(iVariable == Variable::PrecipAcc) {
      return hasVariableCore(Variable::Precip);
   }
   else if(iVariable == Variable::W) {
      return (hasVariableCore(Variable::V) && hasVariableCore(Variable::U)) ||
             (hasVariableCore(Variable::Xwind) && hasVariableCore(Variable::Ywind));
   }
   else if(iVariable == Variable::WD) {
      return (hasVariableCore(Variable::V) && hasVariableCore(Variable::U)) ||
             (hasVariableCore(Variable::Xwind) && hasVariableCore(Variable::Ywind));
   }
   
   // Check if field has been initialized
   std::map<Variable::Type, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
   return it != mFields.end();
}
void File::clear() {
   mFields.clear();
}

long File::getCacheSize() const {
   long size = 0;
   std::map<Variable::Type, std::vector<FieldPtr> >::const_iterator it;
   for(it = mFields.begin(); it != mFields.end(); it++) {
      size += it->second.size() * getNumLat()*getNumLon()*getNumEns()*sizeof(float);
   }
   return size;
}

Uuid File::getUniqueTag() const {
   return mTag;
}
bool File::setLats(vec2 iLats) {
   if(iLats.size() != mNLat || iLats[0].size() != mNLon)
      return false;
   if(mLats != iLats)
      createNewTag();
   mLats = iLats;
   return true;
}
bool File::setLons(vec2 iLons) {
   if(iLons.size() != mNLat || iLons[0].size() != mNLon)
      return false;
   if(mLons != iLons)
      createNewTag();
   mLons = iLons;
   for(int i = 0; i < mLons.size(); i++) {
      for(int j = 0; j < mLons[i].size(); j++) {
         float lon = mLons[i][j];
         if(Util::isValid(lon)) {
            // Ensure lon is between -180 and 180
            int sign = lon / fabs(lon);
            lon = fabs(lon);
            lon = fmod(lon,360); // lon is between 0 and 360
            lon = sign * lon; // lon is between -360 and 306
            if(lon > 180)
               lon = lon - 360;
            else if(lon < -180)
               lon = lon + 360;
            mLons[i][j] = lon;
            assert(mLons[i][j] >= -180.0001 && mLons[i][j] <= 180.0001);
         }
      }
   }

   return true;
}
bool File::setElevs(vec2 iElevs) {
   if(iElevs.size() != mNLat || iElevs[0].size() != mNLon)
      return false;
   mElevs = iElevs;
   return true;
}
bool File::setLandFractions(vec2 iLandFractions) {
   if(iLandFractions.size() != mNLat || iLandFractions[0].size() != mNLon)
      return false;
   mLandFractions = iLandFractions;
   return true;
}
vec2 File::getLats() const {
   return mLats;
}
vec2 File::getLons() const {
   return mLons;
}
vec2 File::getElevs() const {
   return mElevs;
}
vec2 File::getLandFractions() const {
   return mLandFractions;
}
int File::getNumLat() const {
   return mNLat;
}
int File::getNumLon() const {
   return mNLon;
}
int File::getNumEns() const {
   return mNEns;
}
int File::getNumTime() const {
   return mNTime;
}
void File::createNewTag() const {
   mTag = mNextTag; //boost::uuids::random_generator()();
   mNextTag++;
}
void File::setReferenceTime(double iTime) {
   mReferenceTime = iTime;
}
double File::getReferenceTime() const {
   return mReferenceTime;
}
void File::setTimes(std::vector<double> iTimes) {
   if(iTimes.size() != getNumTime()) {
      std::stringstream ss;
      ss << "Setting times array in '" << getFilename() << "' with " << iTimes.size()
         << " elements when the time dimension is " << getNumTime();
      Util::warning(ss.str());
   }
   mTimes = iTimes;
}
std::vector<double> File::getTimes() const {
   return mTimes;
}

std::string File::getDescriptions() {
   std::stringstream ss;
   ss << FileArome::description();
   ss << FileEc::description();
   ss << FilePoint::description();
   ss << FileNorcomQnh::description();
   ss << FileText::description();
   return ss.str();
}
