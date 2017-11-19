#include "Arome.h"
#include <math.h>
#include <netcdf.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileArome::FileArome(std::string iFilename, const Options& iOptions, bool iReadOnly) : FileNetcdfBase(iFilename, iOptions, iReadOnly),
      mXName(""),
      mYName(""),
      mLatName("latitude"),
      mLonName("longitude") {

   iOptions.getValue("lat", mLatName);
   iOptions.getValue("lon", mLonName);

   if(!iOptions.getValue("x", mXName)) {
      if(hasDim("x"))
         mXName = "x";
      else if(hasDim("longitude"))
         mXName = "longitude";
      else if(hasDim("rlon"))
         mXName = "rlon";
      else {
         Util::error("Cannot determine x dimension");
      }
   }
   if(!iOptions.getValue("y", mYName)) {
      if(hasDim("y"))
         mYName = "y";
      else if(hasDim("latitude"))
         mYName = "latitude";
      else if(hasDim("rlat"))
         mYName = "rlat";
      else {
         Util::error("Cannot determine y dimension");
      }
   }

   // Set dimensions
   mNTime = getDimSize("time");
   mNLat  = getDimSize(mYName);
   mNLon  = getDimSize(mXName);
   mNEns  = 1;

   mLats = getLatLonVariable(mLatName);
   mLons = getLatLonVariable(mLonName);
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
      Util::warning("No altitude field available in " + getFilename());
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
   startDataMode();
   // Not cached, retrieve data
   int var = getVar(iVariable);
   int nLat  = mNLat;
   int nLon  = mNLon;

   int numDims = getNumDims(var);

   size_t* count;
   size_t* start;
   long totalCount = nLat*nLon;
   if(numDims == 4) {
      int level = getLevel(iVariable);

      // Variable has a surface dimension
      count = new size_t[4];
      count[0] = 1;
      count[1] = 1;
      count[2] = nLat;
      count[3] = nLon;
      start = new size_t[4];
      start[0] = iTime;
      start[1] = level;
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
   // General rule: Try to define all variables first, then write variables
   // to avoid swapping between define and data mode unnecessarily
   startDefineMode();

   // check if altitudes are valid
   bool isAltitudeValid = false;
   vec2 elevs = getElevs();
   for(int i = 0; i < elevs.size(); i++) {
      for(int j = 0; j < elevs[i].size(); j++) {
         isAltitudeValid = isAltitudeValid || Util::isValid(elevs[i][j]);
      }
   }

   // Define lat/lon/elev
   defineLatLonVariable(mLatName);
   defineLatLonVariable(mLonName);
   if(isAltitudeValid)
      defineLatLonVariable("altitude");

   // Define variables (< 0.01 s)
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      if(!hasVariableCore(varType)) {
         // Create variable
         int dTime    = getDim("time");
         int dLon     = getDim(mXName);
         int dLat     = getDim(mYName);
         int dims[3]  = {dTime, dLat, dLon};
         int var = Util::MV;
         int status = nc_def_var(mFile,variable.c_str(), NC_FLOAT, 3, dims, &var);
         handleNetcdfError(status, "could not define variable");
      }
      int var = getVar(variable);
      float MV = getMissingValue(var); // The output file's missing value indicator
      std::stringstream ss;
      ss << mLonName << " " << mLatName;
      setAttribute(var, "coordinates", ss.str());
      setAttribute(var, "units", Variable::getUnits(varType));
      setAttribute(var, "standard_name", Variable::getStandardName(varType));
      setMissingValue(var, MV);
   }
   defineTimes();
   defineReferenceTime();
   defineGlobalAttributes();

   startDataMode(); // 0.5 seconds
   writeReferenceTime();
   writeTimes(); // 2-3 seconds
   if(isAltitudeValid)
      writeLatLonVariable("altitude");
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      assert(hasVariableCore(varType));
      int var = getVar(variable);
      float MV = getMissingValue(var); // The output file's missing value indicator
      float* values = new float[mNTime*mNLat*mNLon];
      for(int k = 0; k < mNTime*mNLat*mNLon; k++) {
         values[k] = Util::MV;
      }

      // Arrange all data into the 'values' array
      for(int t = 0; t < mNTime; t++) {
         float offset = getOffset(var);
         float scale = getScale(var);
         FieldPtr field = getField(varType, t);
         if(field != NULL) { // TODO: Can't be null if coming from reference
            #pragma omp parallel for
            for(int lat = 0; lat < mNLat; lat++) {
               for(int lon = 0; lon < mNLon; lon++) {
                  float value = (*field)(lat,lon,0);
                  int index = t * mNLat * mNLon + lat * mNLon + lon;
                  if(!Util::isValid(value)) {
                     // Field has missing value indicator and the value is missing
                     // Save values using the file's missing indicator value
                     value = MV;
                  }
                  else {
                     value = ((*field)(lat,lon,0) - offset)/scale;
                  }
                  values[index] = value;
               }
            }
         }
      }

      // Write to file
      int numDims;
      int status = nc_inq_varndims(mFile, var, &numDims);
      handleNetcdfError(status, "could not determine number of dimensions for variable " + variable);
      if(numDims == 4) {
         size_t count[4] = {mNTime, 1, mNLat, mNLon};
         size_t start[4] = {0, 0, 0, 0};
         int status = nc_put_vara_float(mFile, var, start, count, values);
         handleNetcdfError(status, "could not write variable " + variable);
      }
      else if(numDims == 3) {
         size_t count[3] = {mNTime, mNLat, mNLon};
         size_t start[3] = {0, 0, 0};
         int status = nc_put_vara_float(mFile, var, start, count, values);
         handleNetcdfError(status, "could not write variable " + variable);
      }
      else {
         std::stringstream ss;
         ss << "Cannot write variable '" << variable << "' to '" << getFilename() << "'";
         Util::error(ss.str());
      }
      delete[] values;
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
   else if(iVariable == Variable::TD) {
      return "dew_point_temperature_2m";
   }
   else if(iVariable == Variable::Tlevel0) {
      return "air_temperature_ml";
   }
   else if(iVariable == Variable::Tlevel1) {
      return "air_temperature_ml";
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
   else if(iVariable == Variable::PrecipRate) {
      return "lwe_precipitation_rate";
   }
   else if(iVariable == Variable::U) {
      return "eastward_wind_10m";
   }
   else if(iVariable == Variable::Xwind) {
      return "x_wind_10m";
   }
   else if(iVariable == Variable::V) {
      return "northward_wind_10m";
   }
   else if(iVariable == Variable::Ywind) {
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
   // Initialize grid
   vec2 grid;
   grid.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      grid[i].resize(getNumLon());
   }

   if(numDims == 1) {
      // Figure out if this vector has the dimension of lat or lon
      int dLon     = getDim(mXName);
      int dLat     = getDim(mYName);
      int dim;
      int status = nc_inq_vardimid(mFile, var, &dim);
      float* values;
      if(dim == dLon)
         values = new float[getNumLon()];
      else if(dim == dLat)
         values = new float[getNumLat()];
      else
         Util::error("Could not load " + iVar + ". Variables does not have lat or lon dimensions");

      status = nc_get_var_float(mFile, var, values);
      handleNetcdfError(status, "could not get data from variable " + iVar);

      for(int i = 0; i < getNumLat(); i++) {
         for(int j = 0; j < getNumLon(); j++) {
            float value = Util::MV;
            if(dim == dLon)
               value = values[j];
            else
               value = values[i];
            grid[i][j] = value;
         }
      }
      delete[] values;
   }
   else {
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
         Util::error("Cannot read " + iVar + " from AROME file. Must have either 1, 2 or 4 dimensions");
      }
      for(int i = 0; i < getNumLat(); i++) {
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
   }
   return grid;
}
void FileArome::defineLatLonVariable(std::string iVar) {
   if(!hasVar(iVar)) {
      // Create variable
      int dLon     = getDim(mXName);
      int dLat     = getDim(mYName);
      int dims[2]  = {dLat, dLon};
      int var = Util::MV;
      int status = nc_def_var(mFile, iVar.c_str(), NC_FLOAT, 2, dims, &var);
      handleNetcdfError(status, "could not define variable");
   }
}

void FileArome::writeLatLonVariable(std::string iVar) {
   if(!hasVar(iVar)) {
      startDefineMode();
      defineLatLonVariable(iVar);
      startDataMode();
   }
   int var = getVar(iVar);
   float MV = getMissingValue(var);
   int numDims;
   int status = nc_inq_varndims(mFile, var, &numDims);
   handleNetcdfError(status, "could not determine number of dimensions for variable " + iVar);

   // Assign values
   float* values = new float[getNumLon()*getNumLat()];
   vec2 grid;
   if(iVar == "latitude")
      grid = getLats();
   else if(iVar == "longitude")
      grid = getLons();
   else if(iVar == "altitude")
      grid = getElevs();
   for(int i = 0; i < getNumLat(); i++) {
      for(int j = 0; j < getNumLon(); j++) {
         int index = i*getNumLon() + j;
         float value = grid[i][j];
         values[index] = value;
      }
   }

   // Write to file
   if(numDims == 2) {
      status = nc_put_var_float(mFile, var, values);
      handleNetcdfError(status, "could not set data to variable " + iVar);
   }
   else if(numDims == 4) {
      size_t count[4] = {1,1,getNumLat(),getNumLon()};
      size_t start[4] = {0,0,0,0};
      int status = nc_put_vara_float(mFile, var, start, count, values);
      handleNetcdfError(status, "could not set data to variable " + iVar);
   }
   else {
      Util::error("Cannot write " + iVar + " from AROME file. Must have either 2 or 4 dimensions");
   }
   delete[] values;

}
int FileArome::getDate() const {
   return mDate;
}

bool FileArome::isValid(std::string iFilename) {
   bool isValid = false;
   int file;
   int status = nc_open(iFilename.c_str(), NC_NOWRITE, &file);
   if(status == NC_NOERR) {
      isValid = 
               (hasDim(file, "x") || hasDim(file, "rlon")) &&
               (hasDim(file, "y") || hasDim(file, "rlat")) &&
               !hasDim(file, "ensemble_member") &&
               hasVar(file, "latitude") && hasVar(file, "longitude");
      nc_close(file);
   }
   return isValid;
}

int FileArome::getLevel(std::string iVariable) const {
   if(iVariable == "Tlevel0") {
      return 1;
   }
   else if(iVariable == "Tlevel1") {
      return 0;
   }
   return 0;

}

std::string FileArome::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=arome", "AROME file") << std::endl;
   ss << Util::formatDescription("   lat=latitude", "Name of the variable representing latitudes") << std::endl;
   ss << Util::formatDescription("   lon=longitude", "Name of the variable representing longitudes") << std::endl;
   ss << Util::formatDescription("   x=undef", "Name of dimension in the x-direction. If unspecified, the name is auto-detected.") << std::endl;
   ss << Util::formatDescription("   y=undef", "Name of dimension in the y-direction. If unspecified, the name is auto-detected.") << std::endl;
   return ss.str();
}
