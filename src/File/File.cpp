#include "File.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <sstream>
#include "../Util.h"

File::File(std::string iFilename) :
      mFilename(iFilename),
      mTag(boost::uuids::random_generator()()) {
      
}

File* File::getScheme(std::string iFilename, bool iReadOnly) {
   File* file;
   // TODO:
   // Autodetect type based on content
   if(FileArome::isValid(iFilename)) {
      file = new FileArome(iFilename, iReadOnly);
   }
   else if(FileEc::isValid(iFilename)) {
      file = new FileEc(iFilename, iReadOnly);
   }
   else {
      Util::warning("Could not find suitable parser for '" + iFilename + "'");
      return new FileFake(3,3,1,1);
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
         else {
            Util::error("Cannot derive wind speed from variables in file");
         }
      }
      else {
         std::string variableType = Variable::getTypeName(iVariable);
         std::cout << variableType << " not available in '" << getFilename() << "'" << std::endl;
         abort();
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
      return hasVariableCore(Variable::V) && hasVariableCore(Variable::U);
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

boost::uuids::uuid File::getUniqueTag() const {
   return mTag;
}
