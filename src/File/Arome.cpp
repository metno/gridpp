#include "Arome.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileArome::FileArome(std::string iFilename, bool iReadOnly) : FileNetcdf(iFilename, iReadOnly) {
   // Set dimensions
   mNTime = getDimSize("time");
   mNLat  = getDimSize("y");
   mNLon  = getDimSize("x");
   mNEns  = 1;

   mLats = getLatLonVariable("latitude");
   mLons = getLatLonVariable("longitude");
   if(hasVar("surface_geopotential")) {
      FieldPtr elevField = getFieldCore("surface_geopotential", 0);
      mElevs.resize(getNumLat());
      for(int i = 0; i < getNumLat(); i++) {
         mElevs[i].resize(getNumLon());
         for(int j = 0; j < getNumLon(); j++) {
            float value = (*elevField)(i,j,0) / 9.81;
            mElevs[i][j] = value;
         }
      }
      std::cout << "Deriving altitude from geopotential height in " << getFilename() << std::endl;
   }
   else if(hasVar("altitude")) {
      mElevs = getLatLonVariable("altitude");
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

   if(hasVar("land_area_fraction")) {
      mLandFractions = getLatLonVariable("land_area_fraction");
   }
   else {
      mLandFractions.resize(getNumLat());
      for(int i = 0; i < getNumLat(); i++) {
         mLandFractions[i].resize(getNumLon());
         for(int j = 0; j < getNumLon(); j++) {
            mLandFractions[i][j] = Util::MV;
         }
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

FieldPtr FileArome::getFieldCore(Variable::Type iVariable, int iTime) const {
   std::string variableName = getVariableName(iVariable);
   return getFieldCore(variableName, iTime);
}
FieldPtr FileArome::getFieldCore(std::string iVariable, int iTime) const {
   // Not cached, retrieve data
   int var = getVar(iVariable);
   int nLat  = mNLat;
   int nLon  = mNLon;

   int numDims = getNumDims(var);

   size_t* count;
   size_t* start;
   long totalCount = nLat*nLon;
   if(numDims == 4) {
      // Variable has a surface dimension
      count = new size_t[4];
      count[0] = 1;
      count[1] = 1;
      count[2] = nLat;
      count[3] = nLon;
      start = new size_t[4];
      start[0] = iTime;
      start[1] = 0;
      start[2] = 0;
      start[3] = 0;
   }
   else if(numDims == 3) {
      count = new size_t[3];
      count[0] = 1;
      count[1] = nLat;
      count[2] = nLon;
      start = new size_t[3];
      start[0] = iTime;
      start[1] = 0;
      start[2] = 0;
   }
   else {
      std::stringstream ss;
      ss << "Cannot read variable '" << iVariable << "' from '" << getFilename() << "'";
      Util::error(ss.str());
   }
   float* values = new float[nLat*nLon];
   nc_get_vara_float(mFile, var, start, count, values);
   float MV = getMissingValue(var);

   float offset = getOffset(var);
   float scale = getScale(var);
   int index = 0;
   FieldPtr field = getEmptyField();
   for(int lat = 0; lat < nLat; lat++) {
      for(int lon = 0; lon < nLon; lon++) {
         float value = values[index];
         assert(index < totalCount);
         if(value == MV) {
            // Field has missing value indicator and the value is missing
            // Save values using our own internal missing value indicator
            value = Util::MV;
         }
         else {
            value = scale*values[index] + offset;
         }
         (*field)(lat,lon,0) = value;
         index++;
      }
   }
   delete[] values;
   delete[] count;
   delete[] start;
   return field;
}

FileArome::~FileArome() {
}

void FileArome::writeCore(std::vector<Variable::Type> iVariables) {
   int status = ncredef(mFile);
   handleNetcdfError(status, "could not put into define mode");

   // Define variables
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      if(!hasVariableCore(varType)) {
         // Create variable
         int dTime    = getDim("time");
         int dLon     = getDim("x");
         int dLat     = getDim("y");
         int dims[3]  = {dTime, dLat, dLon};
         int var = Util::MV;
         int status = nc_def_var(mFile,variable.c_str(), NC_FLOAT, 3, dims, &var);
         handleNetcdfError(status, "could not define variable");
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
            float* values = new float[mNLat*mNLon];

            int index = 0;
            for(int lat = 0; lat < mNLat; lat++) {
               for(int lon = 0; lon < mNLon; lon++) {
                  float value = (*field)(lat,lon,0);
                  if(!Util::isValid(value)) {
                     // Field has missing value indicator and the value is missing
                     // Save values using the file's missing indicator value
                     value = MV;
                  }
                  else {
                     value = ((*field)(lat,lon,0) - offset)/scale;
                  }
                  values[index] = value;
                  index++;
               }
            }
            int numDims;
            int status = nc_inq_varndims(mFile, var, &numDims);
            handleNetcdfError(status, "could not determine number of dimensions for variable " + variable);
            if(numDims == 4) {
               size_t count[4] = {1, 1, mNLat, mNLon};
               size_t start[4] = {t, 0, 0, 0};
               int status = nc_put_vara_float(mFile, var, start, count, values);
               handleNetcdfError(status, "could not write variable " + variable);
            }
            else if(numDims == 3) {
               size_t count[3] = {1, mNLat, mNLon};
               size_t start[3] = {t, 0, 0};
               int status = nc_put_vara_float(mFile, var, start, count, values);
               handleNetcdfError(status, "could not write variable " + variable);
            }
            else {
               std::stringstream ss;
               ss << "Cannot write variable '" << variable << "' to '" << getFilename() << "'";
               Util::error(ss.str());
            }
            setAttribute(var, "coordinates", "longitude latitude");
            setAttribute(var, "units", Variable::getUnits(varType));
            setAttribute(var, "standard_name", Variable::getStandardName(varType));
            delete[] values;
         }
      }
      setMissingValue(var, MV);
   }
}


std::string FileArome::getVariableName(Variable::Type iVariable) const {
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
   else if(iVariable == Variable::Pop) {
      return "precipitation_amount_prob_low";
   }
   else if(iVariable == Variable::Pop6h) {
      return "precipitation_amount_prob_low_6h";
   }
   else if(iVariable == Variable::PrecipLow) {
      return "precipitation_amount_low_estimate";
   }
   else if(iVariable == Variable::PrecipMiddle) {
      return "precipitation_amount_middle_estimate";
   }
   else if(iVariable == Variable::PrecipHigh) {
      return "precipitation_amount_high_estimate";
   }
   else if(iVariable == Variable::U) {
      return "x_wind_10m";
   }
   else if(iVariable == Variable::V) {
      return "y_wind_10m";
   }
   else if(iVariable == Variable::W) {
      // TODO: Correct name?
      return "windspeed_10m";
   }
   else if(iVariable == Variable::WD) {
      // TODO: Correct name?
      return "winddirection_10m";
   }
   else if(iVariable == Variable::RH) {
      return "relative_humidity_2m";
   }
   else if(iVariable == Variable::Phase) {
      // TODO: Correct name?
      return "phase";
   }
   else if(iVariable == Variable::P) {
      return "surface_air_pressure";
   }
   else if(iVariable == Variable::MSLP) {
      return "air_pressure_at_sea_level";
   }
   else if(iVariable == Variable::QNH) {
      // TODO: What name to use?
      return "qnh";
   }
   else if(iVariable == Variable::SwinAcc) {
      return "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time";
   }
   else if(iVariable == Variable::LwinAcc) {
      return "integral_of_surface_downwelling_longwave_flux_in_air_wrt_time";
   }
   else if(iVariable == Variable::Fake) {
      return "fake";
   }
   else {
      // TODO:
      Util::error("Variable '" + Variable::getTypeName(iVariable) + "' not defined for FileArome");
   }
   return "";
}
vec2 FileArome::getLatLonVariable(std::string iVar) const {
   int var = getVar(iVar);
   float MV = getMissingValue(var);
   int numDims;
   int status = nc_inq_varndims(mFile, var, &numDims);
   handleNetcdfError(status, "could not determine number of dimensions for variable " + iVar);
   float* values = new float[getNumLon()*getNumLat()];
   if(numDims == 2) {
      status = nc_get_var_float(mFile, var, values);
      handleNetcdfError(status, "could not get data from variable " + iVar);
   }
   else if(numDims == 4) {
      size_t count[4] = {1,1,getNumLat(),getNumLon()};
      size_t start[4] = {0,0,0,0};
      int status = nc_get_vara_float(mFile, var, start, count, values);
      handleNetcdfError(status, "could not get data from variable " + iVar);
   }
   else {
      Util::error("Cannot read " + iVar + " from AROME file. Must have either 2 or 4 dimensions");
   }
   vec2 grid;
   grid.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      grid[i].resize(getNumLon());
      for(int j = 0; j < getNumLon(); j++) {
         int index = i*getNumLon() + j;
         float value = values[index];
         if(values[index] == MV)
            value = Util::MV;
         grid[i][j] = value;
         assert(index < getNumLon()*getNumLat());
      }
   }
   delete[] values;
   return grid;
}
int FileArome::getDate() const {
   return mDate;
}

bool FileArome::isValid(std::string iFilename) {
   bool isValid = false;
   int file;
   int status = nc_open(iFilename.c_str(), NC_NOWRITE, &file);
   if(status == NC_NOERR) {
      isValid = hasDim(file, "time") && hasDim(file, "x") && hasDim(file, "y") &&
               !hasDim(file, "ensemble_member") &&
               hasVar(file, "latitude") && hasVar(file, "longitude");
   }
   nc_close(file);
   return isValid;
}

std::string FileArome::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=arome", "AROME file") << std::endl;
   return ss.str();
}
