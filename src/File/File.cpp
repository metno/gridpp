#include "File.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <sstream>
#include "../Util.h"

File::File(std::string iFilename) : mFilename(iFilename) {}

File* File::getScheme(std::string iFilename) {
   File* file;
   // TODO:
   // Autodetect type based on content
   if(FileArome::isValid(iFilename)) {
      file = new FileArome(iFilename);
   }
   else if(FileEc::isValid(iFilename)) {
      file = new FileEc(iFilename);
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
      needsReading = mFields[iVariable][iTime] == NULL;
   }
   else {
      mFields[iVariable].resize(getNumTime());
   }

   if(needsReading) {
      // Load non-derived variable from file
      if(hasVariable(iVariable)) {
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
                     float a1 = (*acc1)[lat][lon][e];
                     float a0 = (*acc0)[lat][lon][e];
                     float value = Util::MV;
                     if(Util::isValid(a1) && Util::isValid(a0)) {
                         value = a1 - a0;
                         if(value < 0)
                            value = 0;
                     }
                     (*field)[lat][lon][e] = value;
                  }
               }
            }
            addField(field, Variable::Precip, t);
         }
      }
      else if(iVariable == Variable::W) {
         if(hasVariable(Variable::U) && hasVariable(Variable::V)) {
            for(int t = 0; t < getNumTime(); t++) {
               FieldPtr windSpeed = getEmptyField();
               const FieldPtr u = getField(Variable::U, t);
               const FieldPtr v = getField(Variable::V, t);
               for(int lat = 0; lat < getNumLat(); lat++) {
                  for(int lon = 0; lon < getNumLon(); lon++) {
                     for(int e = 0; e < getNumEns(); e++) {
                        float currU = (*u)[lat][lon][e];
                        float currV = (*v)[lat][lon][e];
                        (*windSpeed)[lat][lon][e] = sqrt(currU*currU + currV*currV);
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
         std::cout << variableType << " not available in file." << std::endl;
         abort();
      }
   }
   FieldPtr field = mFields[iVariable][iTime];
   return field;
}

File::~File() {
}

void File::write(std::vector<Variable::Type> iVariables) {
   writeCore(iVariables);
}


FieldPtr File::getEmptyField() const {
   return getEmptyField(getNumLat(), getNumLon(), getNumEns());
}
FieldPtr File::getEmptyField(int nLat, int nLon, int nEns) const {
   FieldPtr field = FieldPtr(new Field());
   field->resize(nLat);
   for(int i = 0; i < nLat; i++) {
      (*field)[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         (*field)[i][j].resize(nEns,Util::MV);
      }
   }
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
