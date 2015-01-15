#include "Netcdf.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileNetcdf::FileNetcdf(std::string iFilename) :
   File(iFilename),
   mFile(getFilename().c_str(), NcFile::Write) {
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

FieldPtr FileNetcdf::getFieldCore(Variable::Type iVariable, int iTime) const {
   std::string variable = getVariableName(iVariable);
   // Not cached, retrieve data
   NcVar* var = getVar(variable);
   int nTime = mNTime;
   int nEns  = mNEns;
   int nLat  = mNLat;
   int nLon  = mNLon;

   long count[5] = {1, 1, nEns, nLat, nLon};
   float* values = new float[nTime*1*nEns*nLat*nLon];
   var->set_cur(iTime, 0, 0, 0, 0);
   var->get(values, count);
   float MV = getMissingValue(var);

   float offset = getOffset(var);
   float scale = getScale(var);
   int index = 0;
   FieldPtr field = getEmptyField();
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
            (*field)[lat][lon][e] = value;
            index++;
         }
      }
   }
   delete[] values;
   return field;
}

FileNetcdf::~FileNetcdf() {
   mFile.close();
}

void FileNetcdf::writeCore(std::vector<Variable::Type> iVariables) {
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
         FieldPtr field = getField(varType, t);
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
}


std::string FileNetcdf::getVariableName(Variable::Type iVariable) const {
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

bool FileNetcdf::hasVariable(Variable::Type iVariable) const {
   NcError q(NcError::silent_nonfatal); 
   std::string variable = getVariableName(iVariable);
   NcVar* var = mFile.get_var(variable.c_str());
   return var != NULL;
}

float FileNetcdf::getScale(NcVar* iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* scaleAtt = iVar->get_att("scale_factor");
   float scale  = 1;
   if(scaleAtt != NULL) {
      scale = scaleAtt->as_float(0);
   }
   return scale;
}
float FileNetcdf::getOffset(NcVar* iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* offsetAtt = iVar->get_att("add_offset");
   float offset = 0;
   if(offsetAtt != NULL) {
      offset = offsetAtt->as_float(0);
   }
   return offset;
}

NcDim* FileNetcdf::getDim(std::string iDim) const {
   NcError q(NcError::silent_nonfatal); 
   NcDim* dim = mFile.get_dim(iDim.c_str());
   if(dim == NULL) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return dim;
}
NcVar* FileNetcdf::getVar(std::string iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcVar* var = mFile.get_var(iVar.c_str());
   if(var == NULL) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return var;
}

float FileNetcdf::getMissingValue(const NcVar* iVar) {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* fillValueAtt = iVar->get_att("_FillValue");
   if(fillValueAtt != NULL)
      return fillValueAtt->as_float(0);
   else
      return ncBad_float;//Util::MV;
}
