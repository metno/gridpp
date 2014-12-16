#include "DataFile.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>

DataFile::DataFile(std::string iFilename) :
   mFilename(iFilename), mFile(iFilename.c_str(), NcFile::Write) {
   if(!mFile.is_valid()) {
      std::cout << "Error: Netcdf file " << iFilename << " not valid" << std::endl;
   }

   // Set dimensions
   NcDim* dTime = mFile.get_dim("time");
   NcDim* dEns  = mFile.get_dim("ensemble_member");
   NcDim* dLon  = mFile.get_dim("longitude");
   NcDim* dLat  = mFile.get_dim("latitude");
   mNTime = dTime->size();
   mNEns  = dEns->size();
   mNLat  = dLat->size();
   mNLon  = dLon->size();

   std::cout << "File '" << iFilename << " 'has dimensions " << mNTime << " " << mNEns << " " << mNLat << " " << mNLon << std::endl;
}

const Field& DataFile::getField(Variable::Type iVariable, int iTime) const {
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
                     float value = acc1[lat][lon][e] - acc0[lat][lon][e];
                     if(value < 0)
                        value = 0;
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
   NcError q(NcError::silent_nonfatal); 
   NcVar* var = mFile.get_var(variable.c_str());
   if(var == NULL) {
      std::cout << mFilename << " does not have '" << variable << "'" << std::endl;
      abort();
   }
   int nTime = mNTime;//dTime->size();
   int nEns  = mNEns; //dEns->size();
   int nLat  = mNLat; //dLat->size();
   int nLon  = mNLon; //dLon->size();

   long count[5] = {nTime, 1, nEns, nLat, nLon};
   float* values = new float[nTime*1*nEns*nLat*nLon];
   var->get(values, count);

   float offset = getOffset(var);
   float scale = getScale(var);
   int index = 0;
   for(int t = 0; t < nTime; t++) {
      Field& field = getEmptyField(nLat, nLon, nEns);
      for(int e = 0; e < nEns; e++) {
         for(int lat = 0; lat < nLat; lat++) {
            for(int lon = 0; lon < nLon; lon++) {
               float value = scale*values[index] + offset;
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
         var = mFile.get_var(variable.c_str());
      }
      else {
         // Create variable
         NcDim* dTime = mFile.get_dim("time");
         NcDim* dSurface = mFile.get_dim("surface");
         NcDim* dEns  = mFile.get_dim("ensemble_member");
         NcDim* dLon  = mFile.get_dim("longitude");
         NcDim* dLat  = mFile.get_dim("latitude");
         var = mFile.add_var(variable.c_str(), ncFloat, dTime, dSurface, dEns, dLat, dLon);
      }
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
                     float value = ((*field)[lat][lon][e] - offset)/scale;
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
