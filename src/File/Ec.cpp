#include "Ec.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileEc::FileEc(std::string iFilename, bool iReadOnly) : FileNetcdf(iFilename, iReadOnly) {
   // Set dimensions
   int dTime = getDim("time");
   int dEns  = getDim("ensemble_member");
   int dLon  = getLonDim();
   int dLat  = getLatDim();
   mNTime = getDimSize("time");
   mNEns  = getDimSize("ensemble_member");
   mNLat  = getDimSize(dLat);
   mNLon  = getDimSize(dLon);

   // Retrieve lat/lon/elev
   int vLat = getLatVar();
   int vLon = getLonVar();
   mLats  = getGridValues(vLat);
   mLons  = getGridValues(vLon);

   if(hasVar("altitude")) {
      int vElev = getVar("altitude");
      mElevs = getGridValues(vElev);
   }
   else {
      mElevs.resize(getNumLat());
      for(int i = 0; i < getNumLat(); i++) {
         mElevs[i].resize(getNumLon());
         for(int j = 0; j < getNumLon(); j++) {
            mElevs[i][j] = Util::MV;
         }
      }
   }

   // TODO: No land fraction info in EC files?
   mLandFractions.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      mLandFractions[i].resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         mLandFractions[i][j] = Util::MV;
      }
   }

   if(hasVar("time")) {
      int vTime = getVar("time");
      double* times = new double[mNTime];
      int status = nc_get_var_double(mFile, vTime, times);
      handleNetcdfError(status, "could not get times");
      setTimes(std::vector<double>(times, times+mNTime));
      delete[] times;
   }
   else {
      std::vector<double> times;
      times.resize(getNumTime(), Util::MV);
      setTimes(times);
   }

   if(hasVar("forecast_reference_time")) {
      int vReferenceTime = getVar("forecast_reference_time");
      double referenceTime = getReferenceTime();
      int status = nc_get_var_double(mFile, vReferenceTime, &referenceTime);
      handleNetcdfError(status, "could not get reference time");
      setReferenceTime(referenceTime);
   }

   Util::status( "File '" + iFilename + " 'has dimensions " + getDimenionString());
}

FieldPtr FileEc::getFieldCore(Variable::Type iVariable, int iTime) const {
   std::string variable = getVariableName(iVariable);
   // Not cached, retrieve data
   int var = getVar(variable);
   int nTime = mNTime;
   int nEns  = mNEns;
   int nLat  = mNLat;
   int nLon  = mNLon;

   size_t count[5] = {1, 1, nEns, nLat, nLon};
   size_t start[5] = {iTime, 0, 0, 0, 0};
   float* values = new float[nTime*1*nEns*nLat*nLon];
   nc_get_vara_float(mFile, var, start, count, values);
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
   int status = ncredef(mFile);
   handleNetcdfError(status, "could not put into define mode");

   // Define variables
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      if(variable == "") {
         Util::error("Cannot write variable '" + variable + "' because there EC output file has no definition for it");
      }
      int var = Util::MV;
      if(!hasVariableCore(varType)) {
         // Create variable
         int dTime    = getDim("time");
         int dSurface = getDim("surface");
         int dEns     = getDim("ensemble_member");
         int dLon = getLonDim();
         int dLat = getLatDim();
         int dims[5] = {dTime, dSurface, dEns, dLat, dLon};
         int status = nc_def_var(mFile, variable.c_str(), NC_FLOAT, 5, dims, &var);
         handleNetcdfError(status, "could not define variable '" + variable + "'");
      }
   }
   status = ncendef(mFile);
   handleNetcdfError(status, "could not put into data mode");

   writeTimes();
   writeReferenceTime();
   writeGlobalAttributes();
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      assert(hasVariableCore(varType));
      int var = getVar(variable);
      float MV = getMissingValue(var); // The output file's missing value indicator
      for(int t = 0; t < mNTime; t++) {
         float offset = getOffset(var);
         float scale = getScale(var);
         FieldPtr field = getField(varType, t);
         if(field != NULL) { // TODO: Can't be null if coming from reference
            size_t start[5] = {t, 0, 0, 0, 0};
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
            int numDims = getNumDims(var);
            if(numDims == 5) {
               size_t count[5] = {1,1,mNEns, mNLat, mNLon};
               int status = nc_put_vara_float(mFile, var, start, count, values);
               handleNetcdfError(status, "could not write variable " + variable);
               setAttribute(var, "coordinates", "longitude latitude");
               setAttribute(var, "units", Variable::getUnits(varType));
               setAttribute(var, "standard_name", Variable::getStandardName(varType));
            }
            else {
               std::stringstream ss;
               ss << "Cannot write " << variable << " to '" << getFilename() <<
                             "' because it does not have 5 dimensions. It has " << numDims << " dimensions.";
               Util::warning(ss.str());
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
   else if(iVariable == Variable::W) {
      return "windspeed_10m";
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
   int file;
   int status = nc_open(iFilename.c_str(), NC_NOWRITE, &file);
   if(status == NC_NOERR) {
      isValid = hasDim(file, "time") &&
               (hasVar(file, "lat") || hasVar(file, "latitude")) &&
               (hasVar(file, "lon") || hasVar(file, "longitude")) &&
               hasDim(file, "ensemble_member") &&
               (hasDim(file, "lat") || hasDim(file, "latitude")  || hasDim(file, "y")) &&
               (hasDim(file, "lon") || hasDim(file, "longitude") || hasDim(file, "x"));
      nc_close(file);
   }
   return isValid;
}
vec2 FileEc::getGridValues(int iVar) const {
   // Initialize values
   vec2 grid;
   grid.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      grid[i].resize(getNumLon(), Util::MV);
   }

   // We have a lat/lon grid, where lat/lons are only provided along the pertinent dimension
   // Values are assumed to be constant across the other dimension.
   int numDims = getNumDims(iVar);
   if(numDims == 1) {
      int dim;
      nc_inq_vardimid(mFile, iVar, &dim);
      long size = getDimSize(dim);
      float* values = new float[size];
      nc_get_var_float(mFile, iVar, values);
      // Latitude variable
      if(dim == getLatDim()) {
         for(int i = 0; i < getNumLat(); i++) {
            for(int j = 0; j < getNumLon(); j++) {
               grid[i][j] = values[i];
            }
         }
      }
      // Longitude variable
      else if(dim == getLonDim()) {
         for(int i = 0; i < getNumLat(); i++) {
            for(int j = 0; j < getNumLon(); j++) {
               grid[i][j] = values[j];
            }
         }
      }
      else {
         std::stringstream ss;
         ss << "Missing lat or lon dimension";
         Util::error(ss.str());
      }
      delete[] values;
   }
   // We have a projected grid, where lat and lons are provided for each grid point
   else {
      int N = getNumDims(iVar);
      long count[N];
      int size = 1;
      int indexLat = Util::MV;
      int indexLon = Util::MV;
      int dims[N];
      nc_inq_vardimid(mFile, iVar, dims);
      for(int i = 0; i < N; i++) {
         if(dims[i] == getLatDim()) {
            count[i] = getNumLat();
            size *= count[i];
            indexLat = i;
         }
         else if(dims[i] == getLonDim()) {
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
         ss << "Missing lat and/or lon dimensions";
         Util::error(ss.str());
      }
      float* values = new float[size];
      nc_get_var_float(mFile, iVar, values);
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

int FileEc::getLatDim() const {
   int dLat;
   if(hasDim("y"))
      dLat = getDim("y");
   else if(hasDim("latitude"))
      dLat = getDim("latitude");
   else
      dLat = getDim("lat");
   return dLat;
}
int FileEc::getLonDim() const {
   int dLon;
   if(hasDim("x"))
      dLon = getDim("x");
   else if(hasDim("longitude"))
      dLon = getDim("longitude");
   else
      dLon = getDim("lon");
   return dLon;
}
int FileEc::getLatVar() const {
   int vLat;
   if(hasVar("latitude"))
      vLat = getVar("latitude");
   else
      vLat = getVar("lat");
   return vLat;
}
int FileEc::getLonVar() const {
   int vLon;
   if(hasVar("longitude"))
      vLon = getVar("longitude");
   else
      vLon = getVar("lon");
   return vLon;
}
std::string FileEc::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=ec", "ECMWF ensemble file") << std::endl;
   return ss.str();
}
