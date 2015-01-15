#include "DataFile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "Util.h"

DataFile::DataFile(std::string iFilename) :
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

const Field& DataFile::getField(Variable::Type iVariable, int iTime) const {
   // Check if field is scheduled to be written first
   std::map<Variable::Type, std::vector<Field*> >::const_iterator itW = mWriteFields.find(iVariable);
   if(itW != mWriteFields.end()) {
      if(mWriteFields[iVariable][iTime] != NULL) {
         if(iTime == 1)
            std::cout << "Returning writable field for " << Variable::getTypeName(iVariable) << std::endl;
         return *mWriteFields[iVariable][iTime];
      }
   }

   // Determine if values have been cached
   std::map<Variable::Type, std::vector<Field*> >::const_iterator it = mReadFields.find(iVariable);
   bool needsReading = it == mReadFields.end();
   if(!needsReading) {
      needsReading = mReadFields[iVariable][iTime] == NULL;
   }
   else {
      mReadFields[iVariable].resize(mNTime);
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
   Field* field = mReadFields[iVariable][iTime];
   return *field;
}

void DataFile::loadFields(Variable::Type iVariable) const {
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
void DataFile::saveField(Field* iField, Variable::Type iVariable, int iTime) const {
   if(mReadFields[iVariable][iTime] != NULL)
      delete mReadFields[iVariable][iTime];
   mReadFields[iVariable][iTime] = iField;
}

DataFile::~DataFile() {
   mFile.close();
}


int DataFile::getDate() const {
   return mDate;
}
int DataFile::getInit() const {
   return mInit;
}

void DataFile::write() {
   std::map<Variable::Type, std::vector<Field*> >::const_iterator it;
   for(it = mWriteFields.begin(); it != mWriteFields.end(); it++) {
      std::string variable = getVariableName(it->first);
      NcVar* var;
      if(hasVariable(it->first)) {
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
         Field* field = it->second[t];
         if(field != NULL) {
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


Field& DataFile::getEmptyField() const {
   return getEmptyField(mNLat, mNLon, mNEns);
}
Field& DataFile::getEmptyField(int nLat, int nLon, int nEns) const {
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

void DataFile::addField(Field& iField, Variable::Type iVariable, int iTime) {
   std::map<Variable::Type, std::vector<Field*> >::const_iterator it = mWriteFields.find(iVariable);
   if(it == mWriteFields.end()) {
      mWriteFields[iVariable].resize(mNTime, NULL);
   }

   mWriteFields[iVariable][iTime] = &iField;
}

std::string DataFile::getVariableName(Variable::Type iVariable) const {
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

bool DataFile::hasVariable(Variable::Type iVariable) const {
   NcError q(NcError::silent_nonfatal); 
   std::string variable = getVariableName(iVariable);
   NcVar* var = mFile.get_var(variable.c_str());
   return var != NULL;
}

float DataFile::getScale(NcVar* iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* scaleAtt = iVar->get_att("scale_factor");
   float scale  = 1;
   if(scaleAtt != NULL) {
      scale = scaleAtt->as_float(0);
   }
   return scale;
}
float DataFile::getOffset(NcVar* iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* offsetAtt = iVar->get_att("add_offset");
   float offset = 0;
   if(offsetAtt != NULL) {
      offset = offsetAtt->as_float(0);
   }
   return offset;
}

bool DataFile::hasSameDimensions(const DataFile& iOther) const {
   if(getNumLat() == iOther.getNumLat()
         && getNumLon() == iOther.getNumLon()
         && getNumEns() == iOther.getNumEns()
         && getNumTime() == iOther.getNumTime())
      return true;
   return false;
}

std::string DataFile::getFilename() const {
   return mFilename;
}

std::string DataFile::getDimenionString() const {
   std::stringstream ss;
   ss << "[" << getNumTime() << " " << getNumEns() << " " << getNumLat() << " " << getNumLon()<< "]";
   return ss.str();
}

NcDim* DataFile::getDim(std::string iDim) const {
   NcError q(NcError::silent_nonfatal); 
   NcDim* dim = mFile.get_dim(iDim.c_str());
   if(dim == NULL) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return dim;
}
NcVar* DataFile::getVar(std::string iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcVar* var = mFile.get_var(iVar.c_str());
   if(var == NULL) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return var;
}

float DataFile::getMissingValue(const NcVar* iVar) {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* fillValueAtt = iVar->get_att("_FillValue");
   if(fillValueAtt != NULL)
      return fillValueAtt->as_float(0);
   else
      return Util::MV;
}
