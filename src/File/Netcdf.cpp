#include "Ec.h"
#include <math.h>
#include <netcdf.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileNetcdf::FileNetcdf(std::string iFilename, const Options& iOptions, bool iReadOnly) : File(iFilename, iOptions),
      mInDataMode(true)
{
   int status = nc_open(getFilename().c_str(), iReadOnly ? NC_NOWRITE: NC_WRITE, &mFile);
   if(status != NC_NOERR) {
      Util::error("Could not open NetCDF file " + getFilename());
   }
   // Get defaults
   std::string latVar, lonVar, timeVar, ensDim, latDim, lonDim, timeDim;
   if(!iOptions.getValue("latVar", latVar))
      mLatVar = getLatVar();
   else
      mLatVar = getVar(latVar);

   if(!iOptions.getValue("lonVar", lonVar))
      mLonVar = getLonVar();
   else
      mLonVar = getVar(lonVar);

   if(!iOptions.getValue("timeVar", timeVar))
      mTimeVar = getTimeVar();
   else
      mLonVar = getVar(lonVar);

   // Dimensions
   if(!iOptions.getValue("ensDim", ensDim))
      mEnsDim = getEnsDim();
   else
      mEnsDim = getDim(ensDim);

   if(!iOptions.getValue("latDim", latDim))
      mLatDim = getLatDim();
   else
      mLatDim = getDim(latDim);

   if(!iOptions.getValue("lonDim", lonDim))
      mLonDim = getLonDim();
   else
      mLonDim = getDim(lonDim);

   if(!iOptions.getValue("timeDim", timeDim))
      mTimeDim = getTimeDim();
   else
      mTimeDim = getDim(timeDim);

   // Determine dimension sizes
   mNEns = 1;
   mNLat = 1;
   mNLon = 1;
   mNTime = 1;
   if(Util::isValid(mEnsDim))
      mNEns = getDimSize(mEnsDim);
   if(Util::isValid(mLatDim))
      mNLat = getDimSize(mLatDim);
   if(Util::isValid(mLonDim))
      mNLon = getDimSize(mLonDim);
   if(Util::isValid(mTimeDim))
      mNTime = getDimSize(mTimeDim);

   // Retrieve lat/lon grid
   mLats = getLatLonVariable(mLatVar);
   mLons = getLatLonVariable(mLonVar);

   // Compute elevations
   std::string elevVar;
   mElevVar = Util::MV;
   float elevScale = 1;
   if(!iOptions.getValue("latVar", elevVar)) {
      if(hasVar("altitude")) {
         mElevVar = getVar("altitude");
      }
      else if(hasVar("surface_geopotential")) {
         mElevVar = getVar("surface_geopotential");
         elevScale = 1.0 / 9.81;
      }
   }
   else
      mElevVar = getVar(elevVar);

   if(Util::isValid(mElevVar)) {
      mElevs = getLatLonVariable(mElevVar);
      for(int i = 0; i < getNumLat(); i++) {
         for(int j = 0; j < getNumLon(); j++) {
            mElevs[i][j] *= elevScale;
         }
      }
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


   // Compute times
   if(hasVar("forecast_reference_time")) {
      int vReferenceTime = getVar("forecast_reference_time");
      double referenceTime = getReferenceTime();
      int status = nc_get_var_double(mFile, vReferenceTime, &referenceTime);
      handleNetcdfError(status, "could not get reference time");
      setReferenceTime(referenceTime);
   }

   if(Util::isValid(mTimeVar)) {
      double* times = new double[mNTime];
      int status = nc_get_var_double(mFile, mTimeVar, times);
      handleNetcdfError(status, "could not get times");
      setTimes(std::vector<double>(times, times+mNTime));
      delete[] times;
   }
   else if(Util::isValid(getReferenceTime())) {
      std::stringstream ss;
      ss << "File does not contain times, using forecast reference time as the output time";
      Util::warning(ss.str());
      std::vector<double> times(getNumTime(), getReferenceTime());
      setTimes(times);
   }
   else {
      std::stringstream ss;
      ss << "Cannot determine forecast times, since neither time variable nor reference time variable is available";
      Util::error(ss.str());
   }

   Util::status( "File '" + iFilename + " 'has dimensions " + getDimenionString());
}

FieldPtr FileNetcdf::getFieldCore(Variable::Type iVariable, int iTime) const {
   std::string variableName = getVariableName(iVariable);
   return getFieldCore(variableName, iTime);
}

FieldPtr FileNetcdf::getFieldCore(std::string iVariable, int iTime) const {
   startDataMode();
   int var = getVar(iVariable);
   std::vector<int> dims = getDims(var);

   // Determine which slice to retrieve
   size_t count[dims.size()];
   size_t start[dims.size()];
   int ensPos = Util::MV;
   int latPos = Util::MV;
   int lonPos = Util::MV;
   for(int d = 0; d < dims.size(); d++) {
      int dim = dims[d];
      count[d] = 1;
      start[d] = 0;
      if(dim == mTimeDim) {
         start[d] = iTime;
      }
      else if(dim == mEnsDim) {
         count[d] = getDimSize(dim);
         ensPos = d;
      }
      else if(dim == mLatDim) {
         count[d] = getDimSize(dim);
         latPos = d;
      }
      else if(dim == mLonDim) {
         count[d] = getDimSize(dim);
         lonPos = d;
      }
   }

   // Initialize vector
   size_t size = 1;
   for(int d = 0; d < dims.size(); d++) {
      size *= count[d];
   }
   float* values = new float[size];
   nc_get_vara_float(mFile, var, start, count, values);

   float MV = getMissingValue(var);
   float offset = getOffset(var);
   float scale = getScale(var);

   FieldPtr field = getEmptyField();
   std::vector<int> indices;
   std::vector<int> countVector(count, count + dims.size());
   for(int i = 0; i < size; i++) {
      getIndices(i, countVector, indices);
      int e = 0;
      int lat = 0;
      int lon = 0;
      if(Util::isValid(ensPos))
         e = indices[ensPos];
      if(Util::isValid(latPos))
         lat = indices[latPos];
      if(Util::isValid(lonPos))
         lon = indices[lonPos];
      float value = values[i];
      if(Util::isValid(MV) && value == MV) {
         // Field has missing value indicator and the value is missing
         // Save values using our own internal missing value indicator
         value = Util::MV;
      }
      else {
         value = scale*value + offset;
      }
      (*field)(lat,lon,e) = value;
   }
   delete[] values;
   return field;
}

void FileNetcdf::writeCore(std::vector<Variable::Type> iVariables) {
   startDefineMode();

   // check if altitudes are valid
   bool isAltitudeValid = false;
   vec2 elevs = getElevs();
   for(int i = 0; i < elevs.size(); i++) {
      for(int j = 0; j < elevs[i].size(); j++) {
         isAltitudeValid = isAltitudeValid || Util::isValid(elevs[i][j]);
      }
   }
   if(isAltitudeValid && !hasVar("altitude")) {
      defineAltitude();
   }

   // Define variables
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      std::string typeName = Variable::getTypeName(varType);

      if(variable == "") {
         Util::error("Cannot write variable '" + typeName + "' because there EC output file has no definition for it");
      }
      if(!hasVariableCore(varType)) {
         // Create variable
         int numDims = Util::isValid(mLatDim) + Util::isValid(mLonDim) + Util::isValid(mEnsDim) + Util::isValid(mTimeDim);
         int dims[numDims];
         int counter = 0;
         if(Util::isValid(mTimeDim)) {
            dims[counter] = mTimeDim;
            counter++;
         }
         if(Util::isValid(mEnsDim)) {
            dims[counter] = mEnsDim;
            counter++;
         }
         if(Util::isValid(mLatDim)) {
            dims[counter] = mLatDim;
            counter++;
         }
         if(Util::isValid(mLonDim)) {
            dims[counter] = mLonDim;
            counter++;
         }

         int var = Util::MV;
         int status = nc_def_var(mFile, variable.c_str(), NC_FLOAT, numDims, dims, &var);
         handleNetcdfError(status, "could not define variable '" + variable + "'");
      }
      int var = getVar(variable);
      float MV = getMissingValue(var); // The output file's missing value indicator
      // TODO: Automatically determine if this should be "lon lat" or "longitude latitude"
      setAttribute(var, "coordinates", "lon lat");
      setAttribute(var, "units", Variable::getUnits(varType));
      setAttribute(var, "standard_name", Variable::getStandardName(varType));
   }
   defineTimes();
   defineReferenceTime();
   defineGlobalAttributes();
   startDataMode();

   writeTimes();
   writeReferenceTime();
   if(isAltitudeValid) {
      writeAltitude();
   }
   for(int v = 0; v < iVariables.size(); v++) {
      Variable::Type varType = iVariables[v];
      std::string variable = getVariableName(varType);
      assert(hasVariableCore(varType));
      int var = getVar(variable);
      float MV = getMissingValue(var); // The output file's missing value indicator
      size_t size = 1*1*mNEns*mNLat*mNLon;
      float* values = new float[size];

      std::vector<int> dims = getDims(var);
      size_t count[dims.size()];
      int ensPos = Util::MV;
      int latPos = Util::MV;
      int lonPos = Util::MV;
      int timePos = Util::MV;
      for(int d = 0; d < dims.size(); d++) {
         int dim = dims[d];
         count[d] = 1;
         if(dim == mTimeDim) {
            timePos = d;
         }
         else if(dim == mEnsDim) {
            count[d] = getDimSize(dim);
            ensPos = d;
         }
         else if(dim == mLatDim) {
            count[d] = getDimSize(dim);
            latPos = d;
         }
         else if(dim == mLonDim) {
            count[d] = getDimSize(dim);
            lonPos = d;
         }
      }

      for(int t = 0; t < mNTime; t++) {
         size_t start[dims.size()];
         for(int d = 0; d < dims.size(); d++)
            start[d] = 0;
         if(Util::isValid(timePos))
            start[timePos] = t;

         float offset = getOffset(var);
         float scale = getScale(var);
         FieldPtr field = getField(varType, t);
         if(field != NULL) { // TODO: Can't be null if coming from reference
            std::vector<int> indices;
            std::vector<int> countVector(count, count + dims.size());
            for(int i = 0; i < size; i++) {
               getIndices(i, countVector, indices);
               int e = 0;
               int lat = 0;
               int lon = 0;
               if(Util::isValid(ensPos))
                  e = indices[ensPos];
               if(Util::isValid(latPos))
                  lat = indices[latPos];
               if(Util::isValid(lonPos))
                  lon = indices[lonPos];
               float value = (*field)(lat, lon, e);
               if(Util::isValid(MV) && !Util::isValid(value)) {
                  // Field has missing value indicator and the value is missing
                  // Save values using the file's missing indicator value
                  value = MV;
               }
               else {
                  value = ((*field)(lat,lon,e) - offset)/scale;
               }
               values[i] = value;
            }

            int status = nc_put_vara_float(mFile, var, start, count, values);
            handleNetcdfError(status, "could not write variable " + variable);
         }
      }
      delete[] values;
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
      return "air_temperature_2m";
   }
   else if(iVariable == Variable::TMin) {
      return "air_temperature_2m_min6h";
   }
   else if(iVariable == Variable::TMax) {
      return "air_temperature_2m_max6h";
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

bool FileNetcdf::isValid(std::string iFilename, const Options& iOptions) {
   bool isValid = true;
   int file;
   int status = nc_open(iFilename.c_str(), NC_NOWRITE, &file);
   if(status == NC_NOERR) {
      // Check dimensions
      // isValid = isValid && (hasDim(file, "time") || iOptions.hasValue("timeDim"));
      // isValid = isValid && (hasDim(file, "lat") || hasDim(file, "latitude") || hasDim(file, "y") || iOptions.hasValue("latDim"));
      // isValid = isValid && (hasDim(file, "lon") || hasDim(file, "longitude") || hasDim(file, "x") || iOptions.hasValue("lonDim"));
      // isValid = isValid && (hasDim(file, "ensemble_member") || iOptions.hasValue("ensDim"));

      // Check variables
      // isValid = isValid && (hasVar(file, "time") || iOptions.hasValue("timeVar"));
      isValid = isValid && (hasVar(file, "lat") || hasVar(file, "latitude") || iOptions.hasValue("latVar"));
      isValid = isValid && (hasVar(file, "lon") || hasVar(file, "longitude") || iOptions.hasValue("lonVar"));

      nc_close(file);
   }
   else
      isValid = false;
   return isValid;
}
vec2 FileNetcdf::getGridValues(int iVar) const {
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
      size_t count[N];
      size_t start[N];
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
            size_t dimsize = 1;
            nc_inq_dimlen(mFile, dims[i], &dimsize);
            if(dimsize > 1){
               std::stringstream ss;
               ss << "Lat/lon/elev has an extra non-singleton dimension (dim " << i << ") of length " << dimsize << ". Using index 0 to extract lat/lon/elev.";
               Util::warning(ss.str());
            }
            count[i] = 1;
         }
         start[i] = 0;
      }
      if(!Util::isValid(indexLat) || !Util::isValid(indexLon)) {
         std::stringstream ss;
         ss << "Missing lat and/or lon dimensions";
         Util::error(ss.str());
      }
      float* values = new float[size];
      nc_get_vara_float(mFile, iVar, start, count, values);
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

int FileNetcdf::getLatDim() const {
   int dim = Util::MV;
   if(hasDim("y"))
      dim = getDim("y");
   else if(hasDim("latitude"))
      dim = getDim("latitude");
   else if(hasDim("lat"))
      dim = getDim("lat");
   return dim;
}
int FileNetcdf::getLonDim() const {
   int dim = Util::MV;
   if(hasDim("x"))
      dim = getDim("x");
   else if(hasDim("longitude"))
      dim = getDim("longitude");
   else if(hasDim("lon"))
      dim = getDim("lon");
   return dim;
}
int FileNetcdf::getTimeDim() const {
   int dim = Util::MV;
   if(hasDim("time"))
      dim = getDim("time");
   return dim;
}
int FileNetcdf::getLatVar() const {
   int var = Util::MV;
   if(hasVar("latitude"))
      var = getVar("latitude");
   else if(hasVar("lat"))
      var = getVar("lat");
   else {
      std::stringstream ss;
      ss << "File '" << getFilename() << " does not have a latitude varable";
      Util::error(ss.str());
   }
   return var;
}
int FileNetcdf::getTimeVar() const {
   int var = Util::MV;
   if(hasVar("time"))
      var = getVar("time");
   return var;
}
int FileNetcdf::getLonVar() const {
   int var = Util::MV;
   if(hasVar("longitude"))
      var = getVar("longitude");
   else if(hasVar("lon"))
      var = getVar("lon");
   else {
      std::stringstream ss;
      ss << "File '" << getFilename() << " does not have a longitude varable";
      Util::error(ss.str());
   }
   return var;
}
std::string FileNetcdf::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=ec", "ECMWF ensemble file") << std::endl;
   return ss.str();
}

void FileNetcdf::defineAltitude() {
   // Determine the order of lat and lon dimensions in altitude field
   int indexLat = Util::MV;
   int indexLon = Util::MV;

   // First check if latitude variable has two dimensions
   int vUseDims = Util::MV;
   if(getNumDims(getLatVar()) >= 2) {
      vUseDims = getLatVar();
   }
   else if(hasVar("surface_geopotential")) {
      vUseDims = getVar("surface_geopotential");
   }
   if(Util::isValid(vUseDims)) {
      int N = getNumDims(vUseDims);
      int dimsLat[N];
      nc_inq_vardimid(mFile, vUseDims, dimsLat);
      for(int i = 0; i < N; i++) {
         if(dimsLat[i] == getLatDim())
            indexLat = i;
         if(dimsLat[i] == getLonDim())
            indexLon = i;
      }
   }
   else {
      Util::warning("Could not determine lat/lon ordering when creating altitude variable. Using [lat, lon]");
      indexLat = 0;
      indexLon = 1;
   }

   int dims[2];
   int dLat = getLatDim();
   int dLon = getLonDim();
   if(indexLat < indexLon) {
      dims[0]  = dLat;
      dims[1]  = dLon;
   }
   else {
      dims[0]  = dLon;
      dims[1]  = dLat;
   }
   int var = Util::MV;
   int status = nc_def_var(mFile, "altitude", NC_FLOAT, 2, dims, &var);
   handleNetcdfError(status, "could not define altitude");
}

void FileNetcdf::writeAltitude() const {
   int vElev = getVar("altitude");
   int numDims = getNumDims(vElev);
   if(numDims == 1) {
      Util::error("Cannot write altitude when the variable only has one dimension");
   }
   vec2 elevs = getElevs();

   int N = getNumDims(vElev);
   size_t count[N];
   size_t start[N];
   int size = 1;
   int indexLat = Util::MV;
   int indexLon = Util::MV;
   int dims[N];
   nc_inq_vardimid(mFile, vElev, dims);
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
         size_t dimsize = 1;
         nc_inq_dimlen(mFile, dims[i], &dimsize);
         if(dimsize > 1){
            std::stringstream ss;
            ss << "Lat/lon/elev has an extra non-singleton dimension (dim " << i << ") of length " << dimsize << ". Using index 0 to write lat/lon/elev.";
            Util::warning(ss.str());
         }
         count[i] = 1;
      }
      start[i] = 0;
   }
   float MV = getMissingValue(vElev);
   float* values = new float[size];
   for(int i = 0; i < getNumLat(); i++) {
      for(int j = 0; j < getNumLon(); j++) {
         float elev = elevs[i][j];
         if(Util::isValid(MV) && !Util::isValid(elev))
            elev = MV;
         // Latitude dimension is ordered first
         if(indexLat < indexLon) {
            values[i*getNumLon() + j] = elev;
         }
         // Longitude dimension is ordered first
         else {
            values[j*getNumLat() + i] = elev;
         }
      }
   }
   bool status = nc_put_vara_float(mFile, vElev, start, count, values);
   handleNetcdfError(status, "could not write altitude");
   delete[] values;
}
int FileNetcdf::getVar(std::string iVar) const {
   int var;
   int status = nc_inq_varid(mFile, iVar.c_str(), &var);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return var;
}

int FileNetcdf::getDimSize(std::string iDim) const {
   if(!hasDim(iDim))
      return 0;
   int dim = getDim(iDim);
   size_t len = Util::MV;
   int status = nc_inq_dimlen(mFile, dim, &len);
   return len;
}

int FileNetcdf::getDimSize(int iDim) const {
   size_t len = Util::MV;
   int status = nc_inq_dimlen(mFile, iDim, &len);
   return len;
}

int FileNetcdf::getNumDims(int iVar) const {
   int len = Util::MV;
   int status = nc_inq_varndims(mFile, iVar, &len);
   return len;
}
int FileNetcdf::getEnsDim() const {
   int dim = Util::MV;
   if(hasDim("ensemble_member"))
      dim = getDim("ensemble_member");
   return dim;
}
bool FileNetcdf::hasVar(int iFile, std::string iVar) {
   int var;
   int status = nc_inq_varid(iFile, iVar.c_str(), &var);
   return status == NC_NOERR;
}
bool FileNetcdf::hasVar(std::string iVar) const {
   return hasVar(mFile, iVar);
}
void FileNetcdf::handleNetcdfError(int status, std::string message) const {
   if(status != NC_NOERR) {
      std::stringstream ss;
      if(message == "") {
         ss << "Netcdf error when reading/writing " << getFilename() << ". Netcdf error code: " << status << ".";
      }
      else {
         ss << "Netcdf error for file " << getFilename() << ": " << message << ". "
            << "Netcdf error code: " << status << ".";
      }
      Util::error(ss.str());
   }
}
void FileNetcdf::startDefineMode() const {
   if(mInDataMode) {
      int status = ncredef(mFile);
      handleNetcdfError(status, "could not put into define mode");
   }
   mInDataMode = false;
}
void FileNetcdf::startDataMode() const {
   if(!mInDataMode) {
      int status = ncendef(mFile);
      handleNetcdfError(status, "could not put into data mode");
   }
   mInDataMode = true;
}
std::vector<int> FileNetcdf::getDims(int iVar) const {
   int ndims = getNumDims(iVar);
   assert(ndims >= 0);
   std::vector<int> dims;
   int ids[ndims];
   int status = nc_inq_vardimid(mFile, iVar, ids);
   handleNetcdfError(status, "Could not get dimension");
   for(int i = 0; i < ndims; i++) {
      dims.push_back(ids[i]);
   }
   return dims;
}
float FileNetcdf::getMissingValue(int iVar) const {
   float fillValue;
   int status = nc_get_att_float(mFile, iVar, "_FillValue", &fillValue);
   if(status != NC_NOERR)
      fillValue  = NC_FILL_FLOAT;
   return fillValue;
}
void FileNetcdf::setMissingValue(int iVar, float iValue) const {
   // TODO: Mysterious errors can occur if the existing file was written
   // using an older HDF5 implementation, and the variable has 8 or more
   // attributes already. For more information, see "Corruption Problem
   // In HDF5 1.8.0 through HDF5 1.8.4" on
   // http://www.hdfgroup.org/HDF5/release/known_problems/index.html
   // Does this still apply when using the C-interface to NetCDF?
   if(iValue != NC_FILL_FLOAT) {
      int status = nc_put_att(mFile, iVar, "_FillValue", NC_FLOAT, 1, &iValue);
      handleNetcdfError(status, "could not set missing value flag");
   }
}

float FileNetcdf::getScale(int iVar) const {
   float scale;
   int status = nc_get_att_float(mFile, iVar, "scale_factor", &scale);
   if(status != NC_NOERR)
      scale  = 1;
   return scale;
}
float FileNetcdf::getOffset(int iVar) const {
   float offset;
   int status = nc_get_att_float(mFile, iVar, "add_offset", &offset);
   if(status != NC_NOERR)
      offset  = 0;
   return offset;
}

void FileNetcdf::getIndices(int i, const std::vector<int>& iCount, std::vector<int>& iIndices) const {
   // The last index changes fastest, the first slowest
   int numDims = iCount.size();
   iIndices.resize(numDims);
   int sizeSoFar = 1;
   for(int k = numDims-1; k >= 0; k--) {
      int index = floor(i / sizeSoFar);
      index     = index % iCount[k];
      iIndices[k] = index;
      sizeSoFar *= iCount[k];
   }
}

void FileNetcdf::setAttribute(std::string iVariable, std::string iName, std::string iValue) {
   int var = getVar(iVariable);
   setAttribute(var, iName, iValue);
}

void FileNetcdf::setAttribute(int iVar, std::string iName, std::string iValue) {
   startDefineMode();
   int status = nc_put_att_text(mFile, iVar,iName.c_str(), iValue.size(), iValue.c_str());
   handleNetcdfError(status, "could not set attribute");
}

void FileNetcdf::setGlobalAttribute(std::string iName, std::string iValue) {
   startDefineMode();
   int status = nc_put_att_text(mFile, NC_GLOBAL, iName.c_str(), iValue.size(), iValue.c_str());
   handleNetcdfError(status, "could not set global attribute");
}

void FileNetcdf::appendGlobalAttribute(std::string iName, std::string iValue) {
   int id;
   int status = nc_inq_attid(mFile, NC_GLOBAL, iName.c_str(), &id);
   if(status == NC_NOERR) {
      size_t len;
      int status = nc_inq_attlen(mFile, NC_GLOBAL, iName.c_str(), &len);
      handleNetcdfError(status, "could not determine global attribute length");
      std::stringstream ss;
      if(len > mMaxAttributeLength) {
         std::stringstream ss0;
         ss0 << "Cannot append to attribute with length greater than " << mMaxAttributeLength << ". Resetting to new value";
         Util::warning(ss0.str());
      }
      else {
         char value[len+1];
         status = nc_get_att_text(mFile, NC_GLOBAL, iName.c_str(), value);
         handleNetcdfError(status, "could not append global attribute");
         value[len] = '\0';
         ss << value << "\n";
      }
      ss << iValue;
      setGlobalAttribute(iName, ss.str());
   }
   else {
      setGlobalAttribute(iName, iValue);
   }
}

void FileNetcdf::prependGlobalAttribute(std::string iName, std::string iValue) {
   int id;
   int status = nc_inq_attid(mFile, NC_GLOBAL, iName.c_str(), &id);
   if(status == NC_NOERR) {
      size_t len;
      int status = nc_inq_attlen(mFile, NC_GLOBAL, iName.c_str(), &len);
      handleNetcdfError(status, "could not determine global attribute length");

      std::stringstream ss;
      ss << iValue;
      if(len > mMaxAttributeLength) {
         std::stringstream ss0;
         ss0 << "Cannot prepend to attribute with length greater than " << mMaxAttributeLength << ". Resetting to new value";
         Util::warning(ss0.str());
      }
      else {
         char value[len+1];
         status = nc_get_att_text(mFile, NC_GLOBAL, iName.c_str(), value);
         handleNetcdfError(status, "could not get attribute when prepending a new value");
         value[len] = '\0';
         ss << "\n" << value;
      }
      setGlobalAttribute(iName, ss.str());
   }
   else {
      setGlobalAttribute(iName, iValue);
   }
}

std::string FileNetcdf::getAttribute(std::string iVariable, std::string iName) {
   int var = getVar(iVariable);
   return getAttribute(var, iName);
}

std::string FileNetcdf::getAttribute(int iVar, std::string iName) {
   std::string ret = "";
   int id;
   int status = nc_inq_attid(mFile, iVar, iName.c_str(), &id);
   if(status == NC_NOERR) {
      size_t len;
      int status = nc_inq_attlen(mFile, iVar, iName.c_str(), &len);
      handleNetcdfError(status, "could not determine attribute length");

      if(len > mMaxAttributeLength) {
         std::stringstream ss;
         ss << "Cannot get attribute with length greater than " << mMaxAttributeLength;
         Util::warning(ss.str());
      }
      else {
         char* value = new char[len+1];
         status = nc_get_att_text(mFile, iVar, iName.c_str(), value);
         handleNetcdfError(status, "could not get attribute");
         value[len] = '\0';
         if(status == NC_NOERR)
            ret = std::string(value);
         delete[] value;
      }
   }
   return ret;
}

std::string FileNetcdf::getGlobalAttribute(std::string iName) {
   std::string ret = "";
   int id;
   int status = nc_inq_attid(mFile, NC_GLOBAL, iName.c_str(), &id);
   if(status == NC_NOERR) {
      size_t len;
      int status = nc_inq_attlen(mFile, NC_GLOBAL, iName.c_str(), &len);
      handleNetcdfError(status, "could not determine global attribute length");

      if(len > mMaxAttributeLength) {
         std::stringstream ss;
         ss << "Cannot get attribute with length greater than " << mMaxAttributeLength;
         Util::warning(ss.str());
      }
      else {
         char* value = new char[len+1];
         status = nc_get_att_text(mFile, NC_GLOBAL, iName.c_str(), value);
         handleNetcdfError(status, "could not get global attribute");
         value[len] = '\0';
         if(status == NC_NOERR)
            ret = std::string(value);
         delete[] value;
      }
   }
   return ret;
}

void FileNetcdf::defineTimes() {
   if(!hasVar("time")) {
      if(Util::isValid(mTimeDim)) {
         int id;
         int status = nc_def_var(mFile, "time", NC_DOUBLE, 1, &mTimeDim, &id);
         handleNetcdfError(status, "creating time variable");
      }
      else {
         int id;
         int status = nc_def_var(mFile, "time", NC_DOUBLE, 0, NULL, &id);
         handleNetcdfError(status, "creating time variable");

      }
   }
   int vTime = getVar("time");
   setAttribute(vTime, "long_name", "time");
   setAttribute(vTime, "standard_name", "time");
   setAttribute(vTime, "units", "seconds since 1970-01-01 00:00:00 +00:00");
}
void FileNetcdf::writeTimes() {
   std::vector<double> times = getTimes();
   if(times.size() != getNumTime()) {
      std::stringstream ss;
      ss << "The times specified for NetCDF file '" << getFilename() << "' has " << times.size()
         << " elements, but the time dimension is " << getNumTime() << ". Putting missing values.";
      Util::warning(ss.str());
      times = std::vector<double>(getNumTime(), Util::MV);
   }

   // Convert missing
   for(int i = 0; i < times.size(); i++) {
      if(!Util::isValid(times[i]))
         times[i] = NC_FILL_FLOAT;
   }
   int vTime = getVar("time");
   double timesArr[getNumTime()];
   for(int t = 0; t < getNumTime(); t++) {
      timesArr[t] = times[t];
   }
   int status = nc_put_var_double(mFile, vTime, timesArr);
   handleNetcdfError(status, "could not write times");
}
void FileNetcdf::defineReferenceTime() {
   if(!hasVar("forecast_reference_time")) {
      int id;
      int status = nc_def_var(mFile, "forecast_reference_time", NC_DOUBLE, 0, NULL, &id);
      handleNetcdfError(status, "writing reference time");
   }
   int vTime = getVar("forecast_reference_time");
   setAttribute(vTime, "standard_name", "forecast_reference_time");
   setAttribute(vTime, "units", "seconds since 1970-01-01 00:00:00 +00:00");
}
void FileNetcdf::writeReferenceTime() {
   int vTime = getVar("forecast_reference_time");
   double referenceTime = getReferenceTime();
   if(!Util::isValid(referenceTime))
      referenceTime = NC_FILL_DOUBLE;
   int status = nc_put_var_double(mFile, vTime, &referenceTime);
   handleNetcdfError(status, "could not write reference time");
}
void FileNetcdf::defineGlobalAttributes() {
   if(getGlobalAttribute("Conventions") != "")
      setGlobalAttribute("Conventions", "CF-1.0");
   std::stringstream ss;
   ss << Util::getCurrentTimeStamp() << ": post-processing by gridpp";
   prependGlobalAttribute("history", ss.str());
}
bool FileNetcdf::hasDim(int iFile, std::string iDim) {
   int dim;
   int status = nc_inq_dimid(iFile, iDim.c_str(), &dim);
   return status == NC_NOERR;
}
bool FileNetcdf::hasDim(std::string iDim) const {
   return hasDim(mFile, iDim);
}
bool FileNetcdf::hasVariableCore(Variable::Type iVariable) const {
   std::string variable = getVariableName(iVariable);
   return hasVariableCore(variable);
}
bool FileNetcdf::hasVariableCore(std::string iVariable) const {
   int var=Util::MV;
   int status = nc_inq_varid(mFile, iVariable.c_str(), &var);
   return status == NC_NOERR;
}
vec2 FileNetcdf::getLatLonVariable(int iVar) const {
   float MV = getMissingValue(iVar);
   int numDims;
   std::vector<int> dims = getDims(iVar);
   // Initialize grid
   vec2 grid;
   grid.resize(getNumLat());
   for(int i = 0; i < getNumLat(); i++) {
      grid[i].resize(getNumLon());
   }

   if(dims.size() == 1) {
      // 1D, try to expand vector onto a gri
      float* values;
      if(dims[0] == mLonDim)
         values = new float[getNumLon()];
      else if(dims[0] == mLatDim)
         values = new float[getNumLat()];
      else
         Util::error("Could not load lat/lon. Variable does not have lat or lon dimensions.");

      bool status = nc_get_var_float(mFile, iVar, values);
      handleNetcdfError(status, "could not get data from lat/lon variable");

      for(int i = 0; i < getNumLat(); i++) {
         for(int j = 0; j < getNumLon(); j++) {
            float value = Util::MV;
            if(dims[0] == mLonDim)
               value = values[j];
            else
               value = values[i];
            grid[i][j] = value;
         }
      }
      delete[] values;
   }
   else {
      int latPos = Util::MV;
      int lonPos = Util::MV;
      size_t count[dims.size()];
      size_t start[dims.size()];
      for(int d = 0; d < dims.size(); d++) {
         count[d] = 1;
         start[d] = 0;
         if(dims[d] == mLatDim) {
            latPos = d;
            count[d] = mNLat;
         }
         else if(dims[d] == mLonDim) {
            lonPos = d;
            count[d] = mNLon;
         }
      }
      if(!Util::isValid(latPos) || !Util::isValid(lonPos)) {
         std::stringstream ss;
         ss << "lat/lon field does not have lat or lon dimension";
         Util::error(ss.str());
      }

      float* values = new float[getNumLon()*getNumLat()];
      bool status = nc_get_vara_float(mFile, iVar, start, count, values);
      handleNetcdfError(status, "could not get data from lat/lon variable");
      for(int i = 0; i < getNumLat(); i++) {
         for(int j = 0; j < getNumLon(); j++) {
            int index = 0;
            if(latPos < lonPos) {
               index = i*getNumLon() + j;
            }
            else {
               index = j*getNumLat() + i;
            }
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
