#include "File.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "Util.h"

File::File(std::string iFilename) :
   mFilename(iFilename), mFile(iFilename.c_str(), NcFile::Write) {
   if(!mFile.is_valid()) {
      Util::error("Error: Netcdf file " + iFilename + " not valid");
   }

   // Set dimensions
   NcDim* dTime = getDim("time");
   NcDim* dEns  = getDim("ensemble_member");
   NcDim* dLon  = getDim("longitude");
   NcDim* dLat  = getDim("latitude");
   mNTime = dTime->size();
   mNEns  = dEns->size();
   mNLat  = dLat->size();
   mNLon  = dLon->size();

   Util::status( "File '" + iFilename + " 'has dimensions " + getDimenionString());
}

Field& File::getField(Variable::Type iVariable, int iTime) const {
   // Determine if values have been cached
   std::map<Variable::Type, std::vector<Field*> >::const_iterator it = mFields.find(iVariable);
   bool needsReading = it == mFields.end();
   if(!needsReading) {
      needsReading = mFields[iVariable][iTime] == NULL;
   }
   else {
      mFields[iVariable].resize(mNTime);
   }

   if(needsReading) {
      // Load non-derived variable from file
      if(hasVariable(iVariable)) {
         loadFields(iVariable);
      }
      // Try to derive the field
      else if(iVariable == Variable::Precip) {
         // Deaccumulate
         Field& field = getEmptyField();
         saveField(&field, Variable::Precip, 0); // First offset is 0

         for(int t = 1; t < mNTime; t++) {
            Field& field = getEmptyField();
            const Field& acc0  = getField(Variable::PrecipAcc, t-1);
            const Field& acc1  = getField(Variable::PrecipAcc, t);
            for(int lat = 0; lat < mNLat; lat++) {
               for(int lon = 0; lon < mNLon; lon++) {
                  for(int e = 0; e < mNEns; e++) {
                     float a1 = acc1[lat][lon][e];
                     float a0 = acc0[lat][lon][e];
                     float value = Util::MV;
                     if(Util::isValid(a1) && Util::isValid(a0)) {
                         value = a1 - a0;
                         if(value < 0)
                            value = 0;
                     }
                     field[lat][lon][e] = value;
                  }
               }
            }
            saveField(&field, Variable::Precip, t);
         }
      }
      else {
         std::string variableType = Variable::getTypeName(iVariable);
         std::cout << "File " << mFilename << " does not contain a recognized "
                   << variableType << " variable" << std::endl;
         abort();
      }
   }
   Field* field = mFields[iVariable][iTime];
   return *field;
}

void File::loadFields(Variable::Type iVariable) const {
   std::string variable = getVariableName(iVariable);
   // Not cached, retrieve data
   NcVar* var = getVar(variable);
   int nTime = mNTime;//dTime->size();
   int nEns  = mNEns; //dEns->size();
   int nLat  = mNLat; //dLat->size();
   int nLon  = mNLon; //dLon->size();

   long count[5] = {nTime, 1, nEns, nLat, nLon};
   float* values = new float[nTime*1*nEns*nLat*nLon];
   var->get(values, count);
   float MV = getMissingValue(var);

   float offset = getOffset(var);
   float scale = getScale(var);
   int index = 0;
   for(int t = 0; t < nTime; t++) {
      Field& field = getEmptyField(nLat, nLon, nEns);
      for(int e = 0; e < nEns; e++) {
         for(int lat = 0; lat < nLat; lat++) {
            for(int lon = 0; lon < nLon; lon++) {
               float value = values[index];
               if(Util::isValid(MV) && value == MV) {
                  // Field has missing value indicator and the value is missing
                  // Save values using our own internal missing value indicator
                  value = Util::MV;
               }
               else {
                  value = scale*values[index] + offset;
               }
               field[lat][lon][e] = value;
               index++;
            }
         }
      }
      saveField(&field, iVariable, t);
   }
   delete[] values;
}
void File::saveField(Field* iField, Variable::Type iVariable, int iTime) const {
   if(mFields[iVariable][iTime] != NULL)
      delete mFields[iVariable][iTime];
   mFields[iVariable][iTime] = iField;
}

File::~File() {
   mFile.close();
}

int File::getDate() const {
   return mDate;
}
int File::getInit() const {
   return mInit;
}

void File::write(std::vector<Variable::Type> iVariables) {
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      NcVar* var;
      if(hasVariable(varType)) {
         var = getVar(variable);
      }
      else {
         // Create variable
         NcDim* dTime    = getDim("time");
         NcDim* dSurface = getDim("surface");
         NcDim* dEns     = getDim("ensemble_member");
         NcDim* dLon     = getDim("longitude");
         NcDim* dLat     = getDim("latitude");
         var = mFile.add_var(variable.c_str(), ncFloat, dTime, dSurface, dEns, dLat, dLon);
      }
      float MV = getMissingValue(var); // The output file's missing value indicator
      for(int t = 0; t < mNTime; t++) {
         float offset = getOffset(var);
         float scale = getScale(var);
         Field* field = &getField(varType, t);
         if(field != NULL) { // TODO: Can't be null if coming from reference
            var->set_cur(t, 0, 0, 0, 0);
            float* values = new float[mNTime*1*mNEns*mNLat*mNLon];

            int index = 0;
            for(int e = 0; e < mNEns; e++) {
               for(int lat = 0; lat < mNLat; lat++) {
                  for(int lon = 0; lon < mNLon; lon++) {
                     float value = (*field)[lat][lon][e];
                     if(Util::isValid(MV) && !Util::isValid(value)) {
                        // Field has missing value indicator and the value is missing
                        // Save values using the file's missing indicator value
                        value = MV;
                     }
                     else {
                        value = ((*field)[lat][lon][e] - offset)/scale;
                     }
                     values[index] = value;
                     index++;
                  }
               }
            }
            var->put(values, 1, 1, mNEns, mNLat, mNLon);
         }
      }
   }
   mFile.close();
}


Field& File::getEmptyField() const {
   return getEmptyField(mNLat, mNLon, mNEns);
}
Field& File::getEmptyField(int nLat, int nLon, int nEns) const {
   Field* field = new Field();
   field->resize(nLat);
   for(int i = 0; i < nLat; i++) {
      (*field)[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         (*field)[i][j].resize(nEns,0);
      }
   }
   return *field;
}

void File::addField(Field& iField, Variable::Type iVariable, int iTime) {
   std::map<Variable::Type, std::vector<Field*> >::const_iterator it = mFields.find(iVariable);
   if(it == mFields.end()) {
      mFields[iVariable].resize(mNTime, NULL);
   }

   mFields[iVariable][iTime] = &iField;
}

std::string File::getVariableName(Variable::Type iVariable) const {
   if(iVariable == Variable::PrecipAcc) {
      return "precipitation_amount_acc";
   }
   else if(iVariable == Variable::Cloud) {
      return "cloud_area_fraction";
   }
   else if(iVariable == Variable::T) {
      return "t";
   }
   else if(iVariable == Variable::Precip) {
      return "precipitation_amount";
   }
   return "";
}

bool File::hasVariable(Variable::Type iVariable) const {
   NcError q(NcError::silent_nonfatal); 
   std::string variable = getVariableName(iVariable);
   NcVar* var = mFile.get_var(variable.c_str());
   return var != NULL;
}

float File::getScale(NcVar* iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* scaleAtt = iVar->get_att("scale_factor");
   float scale  = 1;
   if(scaleAtt != NULL) {
      scale = scaleAtt->as_float(0);
   }
   return scale;
}
float File::getOffset(NcVar* iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* offsetAtt = iVar->get_att("add_offset");
   float offset = 0;
   if(offsetAtt != NULL) {
      offset = offsetAtt->as_float(0);
   }
   return offset;
}

bool File::hasSameDimensions(const File& iOther) const {
   if(getNumLat() == iOther.getNumLat()
         && getNumLon() == iOther.getNumLon()
         && getNumEns() == iOther.getNumEns()
         && getNumTime() == iOther.getNumTime())
      return true;
   return false;
}

std::string File::getFilename() const {
   return mFilename;
}

std::string File::getDimenionString() const {
   std::stringstream ss;
   ss << "[" << getNumTime() << " " << getNumEns() << " " << getNumLat() << " " << getNumLon()<< "]";
   return ss.str();
}

NcDim* File::getDim(std::string iDim) const {
   NcError q(NcError::silent_nonfatal); 
   NcDim* dim = mFile.get_dim(iDim.c_str());
   if(dim == NULL) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return dim;
}
NcVar* File::getVar(std::string iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcVar* var = mFile.get_var(iVar.c_str());
   if(var == NULL) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return var;
}

float File::getMissingValue(const NcVar* iVar) {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* fillValueAtt = iVar->get_att("_FillValue");
   if(fillValueAtt != NULL)
      return fillValueAtt->as_float(0);
   else
      return Util::MV;
}
