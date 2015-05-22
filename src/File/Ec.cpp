#include "Ec.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileEc::FileEc(std::string iFilename, bool iReadOnly) : FileNetcdf(iFilename, iReadOnly) {
   // Set dimensions
   NcDim* dTime = getDim("time");
   NcDim* dEns  = getDim("ensemble_member");
   NcDim* dLon  = getLonDim();
   NcDim* dLat  = getLatDim();
   mNTime = dTime->size();
   mNEns  = dEns->size();
   mNLat  = dLat->size();
   mNLon  = dLon->size();

   // Retrieve lat/lon/elev
   NcVar* vLat = getLatVar();
   NcVar* vLon = getLonVar();
   NcVar* vElev = getVar("altitude");
   mLats  = getGridValues(vLat);
   mLons  = getGridValues(vLon);
   mElevs = getGridValues(vElev);

   if(hasVar("time")) {
      NcVar* vTime = getVar("time");
      double* times = new double[mNTime];
      vTime->get(times , mNTime);
      setTimes(std::vector<double>(times, times+mNTime));
      delete[] times;
   }
   else {
      std::vector<double> times;
      times.resize(getNumTime(), Util::MV);
      setTimes(times);
   }

   if(hasVar("forecast_reference_time")) {
      NcVar* vReferenceTime = getVar("forecast_reference_time");
      double referenceTime = getReferenceTime();
      vReferenceTime->get(&referenceTime, 1);
      setReferenceTime(referenceTime);
   }

   Util::status( "File '" + iFilename + " 'has dimensions " + getDimenionString());
}

FieldPtr FileEc::getFieldCore(Variable::Type iVariable, int iTime) const {
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
            (*field)(lat,lon,e) = value;
            index++;
         }
      }
   }
   delete[] values;
   return field;
}

void FileEc::writeCore(std::vector<Variable::Type> iVariables) {
   writeTimes();
   writeReferenceTime();
   writeGlobalAttributes();
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      NcVar* var;
      if(hasVariableCore(varType)) {
         var = getVar(variable);
      }
      else {
         // Create variable
         NcDim* dTime    = getDim("time");
         NcDim* dSurface = getDim("surface");
         NcDim* dEns     = getDim("ensemble_member");
         NcDim* dLon     = getLonDim();
         NcDim* dLat     = getLatDim();
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
                     float value = (*field)(lat,lon,e);
                     if(Util::isValid(MV) && !Util::isValid(value)) {
                        // Field has missing value indicator and the value is missing
                        // Save values using the file's missing indicator value
                        value = MV;
                     }
                     else {
                        value = ((*field)(lat,lon,e) - offset)/scale;
                     }
                     values[index] = value;
                     index++;
                  }
               }
            }
            if(var->num_dims() == 5) {
               var->put(values, 1, 1, mNEns, mNLat, mNLon);
               setAttribute(var, "coordinates", "longitude latitude");
               setAttribute(var, "units", Variable::getUnits(varType));
               setAttribute(var, "standard_name", Variable::getStandardName(varType));
            }
            else {
               Util::warning("Cannot write " + variable + " to '" + getFilename() +
                             "' because it does not have 5 dimensions");
            }
         }
      }
   }
}


std::string FileEc::getVariableName(Variable::Type iVariable) const {
   if(iVariable == Variable::PrecipAcc) {
      return "precipitation_amount_acc";
   }
   else if(iVariable == Variable::Cloud) {
      return "cloud_area_fraction";
   }
   else if(iVariable == Variable::T) {
      return "air_temperature_2m";
   }
   else if(iVariable == Variable::Precip) {
      return "precipitation_amount";
   }
   else if(iVariable == Variable::U) {
      return "x_wind_10m";
   }
   else if(iVariable == Variable::V) {
      return "y_wind_10m";
   }
   else if(iVariable == Variable::MSLP) {
      return "sea_level_pressure";
   }
   else if(iVariable == Variable::RH) {
      return "relative_humidity_2m";
   }
   return "";
}

bool FileEc::isValid(std::string iFilename) {
   bool isValid = false;
   NcFile file = NcFile(iFilename.c_str(), NcFile::ReadOnly);
   if(file.is_valid()) {
      isValid = hasDim(file, "time") &&
               (hasVar(file, "lat") || hasVar(file, "latitude")) &&
               (hasVar(file, "lon") || hasVar(file, "longitude")) &&
               hasDim(file, "ensemble_member") &&
               (hasDim(file, "lat") || hasDim(file, "latitude")  || hasDim(file, "y")) &&
               (hasDim(file, "lon") || hasDim(file, "longitude") || hasDim(file, "x"));
      file.close();
   }
   return isValid;
}
vec2 FileEc::getGridValues(NcVar* iVar) const {
   // Initialize values
   vec2 grid;
   grid.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      grid[i].resize(getNumLon(), Util::MV);
   }

   // We have a lat/lon grid, where lat/lons are only provided along the pertinent dimension
   // Values are assumed to be constant across the other dimension.
   if(iVar->num_dims() == 1) {
      long size = iVar->get_dim(0)->size();
      long count[1] = {size};
      float* values = new float[size];
      iVar->get(values, count);
      // Latitude variable
      if(iVar->get_dim(0) == getLatDim()) {
         for(int i = 0; i < getNumLat(); i++) {
            for(int j = 0; j < getNumLon(); j++) {
               grid[i][j] = values[i];
            }
         }
      }
      // Longitude variable
      else if(iVar->get_dim(0) == getLonDim()) {
         for(int i = 0; i < getNumLat(); i++) {
            for(int j = 0; j < getNumLon(); j++) {
               grid[i][j] = values[j];
            }
         }
      }
      else {
         std::stringstream ss;
         ss << "Variable " << iVar->name() << " does not have lat or lon dimension";
         Util::error(ss.str());
      }
      delete[] values;
   }
   // We have a projected grid, where lat and lons are provided for each grid point
   else {
      int N = iVar->num_dims();
      long count[N];
      int size = 1;
      int indexLat = Util::MV;
      int indexLon = Util::MV;
      for(int i = 0; i < N; i++) {
         if(iVar->get_dim(i) == getLatDim()) {
            count[i] = getNumLat();
            size *= count[i];
            indexLat = i;
         }
         else if(iVar->get_dim(i) == getLonDim()) {
            count[i] = getNumLon();
            size *= count[i];
            indexLon = i;
         }
         else {
            count[i] = 1;
         }
      }
      if(!Util::isValid(indexLat) || !Util::isValid(indexLon)) {
         std::stringstream ss;
         ss << "Variable " << iVar->name() << " does not have lat and/or lon dimensions";
         Util::error(ss.str());
      }
      float* values = new float[size];
      iVar->get(values, count);
      for(int i = 0; i < getNumLat(); i++) {
         for(int j = 0; j < getNumLon(); j++) {
            // Latitude dimension is ordered first
            if(indexLat < indexLon) {
               grid[i][j] = values[i*getNumLon() + j];
            }
            // Longitude dimension is ordered first
            else {
               grid[i][j] = values[j*getNumLat() + i];
            }
         }
      }
      delete[] values;
   }
   return grid;
}

NcDim* FileEc::getLatDim() const {
   NcDim* dLat;
   if(hasDim("y"))
      dLat = getDim("y");
   else if(hasDim("latitude"))
      dLat = getDim("latitude");
   else
      dLat = getDim("lat");
   return dLat;
}
NcDim* FileEc::getLonDim() const {
   NcDim* dLon;
   if(hasDim("x"))
      dLon = getDim("x");
   else if(hasDim("longitude"))
      dLon = getDim("longitude");
   else
      dLon = getDim("lon");
   return dLon;
}
NcVar* FileEc::getLatVar() const {
   NcVar* vLat;
   if(hasVar("latitude"))
      vLat = getVar("latitude");
   else
      vLat = getVar("lat");
   return vLat;
}
NcVar* FileEc::getLonVar() const {
   NcVar* vLon;
   if(hasVar("longitude"))
      vLon = getVar("longitude");
   else
      vLon = getVar("lon");
   return vLon;
}
