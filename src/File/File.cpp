#include "File.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include "../Util.h"
#include "../Options.h"
Uuid File::mNextTag = 0;

File::File(std::string iFilename, const Options& iOptions) :
      mFilename(iFilename),
      mReferenceTime(Util::MV) {
   createNewTag();

   std::string variablesFile;
   if(!iOptions.getValue("variables", variablesFile))
      mVariables = VariableMap::loadDefaults();
   else
      mVariables = VariableMap::load(variablesFile);
}

File* File::getScheme(std::string iFilename, const Options& iOptions, bool iReadOnly) {
   File* file;
   // Determine the filetype, either through user-specified option type=...
   // or by autodetecting.
   std::string type = "";
   if(!iOptions.getValue("type", type)) {
      // Autodetect type based on content
      if(FileNetcdf::isValid(iFilename, iOptions)) {
         type = "netcdf";
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
   else if(type == "netcdf") {
      file = new FileNetcdf(iFilename, iOptions, iReadOnly);
   }
   /*
    * TODO:
   else if(type == "point") {
      file = new FilePoint(iFilename, iOptions);
   }
   else if(type == "text") {
      file = new FileText(iFilename, iOptions);
   }
   else if(type == "norcomQnh") {
      file = new FileNorcomQnh(iFilename, iOptions);
   }
  */
   else {
      Util::error("Could not understand file type " + type);
   }
   return file;
}

FieldPtr File::getField(Variable::Type iVariable, int iTime) const {
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].getType() == iVariable) {
         return getField(mVariables[i], iTime);
      }
   }
   // TODO:
   Util::error("Could not get field");
}

FieldPtr File::getField(std::string iVariable, int iTime) const {
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].getName() == iVariable) {
         return getField(mVariables[i], iTime);
      }
   }
   // TODO:
   Util::error("Could not get field");
}

FieldPtr File::getField(const Variable& iVariable, int iTime) const {
   // Determine if values have been cached
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
   bool needsReading = it == mFields.end();
   if(!needsReading) {
      if(mFields[iVariable].size() <= iTime) {
         std::stringstream ss;
         ss << "Attempted to access variable '" << iVariable.getName() << "' for time " << iTime
            << " in file '" << getFilename() << "'";
         Util::error(ss.str());
      }

      needsReading = mFields[iVariable][iTime] == NULL;
   }
   else {
      if(getNumTime() <= iTime) {
         std::stringstream ss;
         ss << "Attempted to access variable '" << iVariable.getName() << "' for time " << iTime
            << " in file '" << getFilename() << "'";
         Util::error(ss.str());
      }
      mFields[iVariable].resize(getNumTime());
   }

   if(needsReading) {
      // Load non-derived variable from file
      if(hasVariableCore(iVariable)) {
         mFields[iVariable][iTime] = getFieldCore(iVariable, iTime);
      }
      // Try to derive the field
      else if(iVariable.getType() == Variable::Precip) {
         Variable outputVariable;
         bool status = getVariable(Variable::PrecipAcc, outputVariable);
         if(!status) {
            Util::error("Cannot derive Precip");
         }
         // Deaccumulate
         FieldPtr field = getEmptyField();
         addField(field, iVariable, 0); // First offset is 0

         for(int t = 1; t < getNumTime(); t++) {
            FieldPtr field = getEmptyField();
            const FieldPtr acc0  = getField(outputVariable, t-1);
            const FieldPtr acc1  = getField(outputVariable, t);
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
            addField(field, iVariable, t);
         }
      }
      else if(iVariable.getType() == Variable::PrecipAcc) {
         Variable outputVariable;
         bool status = getVariable(Variable::Precip, outputVariable);
         if(!status) {
            Util::error("Cannot derive PrecipAcc");
         }
         // Accumulate
         FieldPtr prevAccum = getEmptyField(0);
         addField(prevAccum, iVariable, 0); // First offset is 0

         for(int t = 1; t < getNumTime(); t++) {
            FieldPtr currAccum = getEmptyField();
            const FieldPtr currPrecip  = getField(outputVariable, t);
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
            addField(currAccum, iVariable, t);
            prevAccum = currAccum;
         }
      }
      else if(iVariable.getType() == Variable::W) {
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
                        if(Util::isValid(currU) && Util::isValid(currV)) {
                           (*windSpeed)(lat,lon,e) = sqrt(currU*currU + currV*currV);
                        }
                     }
                  }
               }
               addField(windSpeed, iVariable, t);
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
                        if(Util::isValid(currX) && Util::isValid(currY)) {
                           (*windSpeed)(lat,lon,e) = sqrt(currX*currX + currY*currY);
                        }
                     }
                  }
               }
               addField(windSpeed, iVariable, t);
            }
         }
         else {
            Util::error("Cannot derive wind speed from variables in file");
         }
      }
      else if(iVariable.getType() == Variable::WD) {
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
                        if(Util::isValid(currU) && Util::isValid(currV)) {
                           float dir = std::atan2(-currU,-currV) * 180 / Util::pi;
                           if(dir < 0)
                              dir += 360;
                           (*windDir)(lat,lon,e) = dir;
                        }
                     }
                  }
               }
               addField(windDir, iVariable, t);
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
                        if(Util::isValid(currX) && Util::isValid(currY)) {
                           float dir = std::atan2(-currX,-currY) * 180 / Util::pi;
                           if(dir < 0)
                              dir += 360;
                           (*windDir)(lat,lon,e) = dir;
                        }
                     }
                  }
               }
               addField(windDir, iVariable, t);
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
         std::string variableType = iVariable.getName();
         Util::warning(variableType + " not available in '" + getFilename() + "'");
         for(int t = 0; t < getNumTime(); t++) {
            FieldPtr field = getEmptyField();
            addField(field, iVariable, t);
         }
      }
   }
   if(mFields[iVariable].size() <= iTime) {
      std::stringstream ss;
      ss << "Attempted to access variable '" << iVariable.getName() << "' for time " << iTime
         << " in file '" << getFilename() << "'";
      Util::error(ss.str());
   }
   FieldPtr field = mFields[iVariable][iTime];
   return field;
}

File::~File() {
}

void File::write(std::vector<Variable> iVariables) {
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

void File::addField(FieldPtr iField, const Variable& iVariable, int iTime) const {
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
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

void File::initNewVariable(const Variable& iVariable) {
   if(!hasVariable(iVariable)) {
      for(int t = 0; t < getNumTime(); t++) {
         addField(getEmptyField(), iVariable, t);
      }
   }
}
void File::setVariables(std::vector<Variable> iVariables) {
   mVariables = iVariables;
}
bool File::hasVariable(std::string iVariable) const {
   // TODO: Check derived variables
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].getName() == iVariable)
         return true;
   }
   return false;

}
bool File::hasVariable(Variable::Type iVariable) const {
   // TODO
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].getType() == iVariable)
         return true;
   }
   return false;
}
bool File::hasVariable(const Variable& iVariable) const {
   bool status = hasVariableCore(iVariable);
   if(status)
      return true;
// Check if field is derivable
   if(iVariable.getType() == Variable::Precip) {
      return hasVariableCore(Variable::PrecipAcc);
   }
   else if(iVariable.getType() == Variable::PrecipAcc) {
      return hasVariableCore(Variable::Precip);
   }
   else if(iVariable.getType() == Variable::W) {
      return (hasVariableCore(Variable::V) && hasVariableCore(Variable::U)) ||
             (hasVariableCore(Variable::Xwind) && hasVariableCore(Variable::Ywind));
   }
   else if(iVariable.getType() == Variable::WD) {
      return (hasVariableCore(Variable::V) && hasVariableCore(Variable::U)) ||
             (hasVariableCore(Variable::Xwind) && hasVariableCore(Variable::Ywind));
   }

   // Check if field has been initialized
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
   return it != mFields.end();
}

bool File::hasVariableWithoutDeriving(Variable::Type iVariable) const {
   bool status = hasVariableCore(iVariable);
   if(status)
      return true;

   // Check if field has been initialized
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it;
   for(it == mFields.begin(); it != mFields.end(); it++) {
      if(it->first.getType() == iVariable)
         return true;
   }
   return false;
}
void File::clear() {
   mFields.clear();
}

long File::getCacheSize() const {
   long size = 0;
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it;
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
      mNTime = iTimes.size();
   }
   mTimes = iTimes;
}
std::vector<double> File::getTimes() const {
   return mTimes;
}

bool File::hasVariableCore(Variable::Type iVariable) const {
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].getType() == iVariable) {
         return true;
      }
   }
   return false;
}

bool File::hasVariableCore(std::string iVariable) const {
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].getName() == iVariable) {
         return true;
      }
   }
   return false;
}

std::string File::getVariableName(const Variable& iVariable) const {
   // TODO:
   //if(mVariableMap.has(iVariable))
   //   return mVariableMap.get(iVariable);
   //else
      return "";
}

bool File::getVariable(Variable::Type iVariableType, Variable& iVariable) const {
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].getType() == iVariableType) {
         iVariable = mVariables[i];
         return true;
      }
   }
   return false;
}

std::string File::getDescriptions() {
   std::stringstream ss;
   ss << FileNetcdf::description();
   ss << FilePoint::description();
   ss << FileNorcomQnh::description();
   ss << FileText::description();
   return ss.str();
}
