#include "Netcdf.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include "../NetcdfUtil.h"
#include <assert.h>
#include <set>
#include <fstream>
#include <math.h>

ParameterFileNetcdf::ParameterFileNetcdf(const Options& iOptions) : ParameterFile(iOptions),
      mDimName("coeff"),
      mVarName("coefficient") {
   int file;
   int status = nc_open(getFilename().c_str(), NC_NOWRITE, &file);
   handleNetcdfError(status, "invalid parameter file");
   iOptions.getValue("dimName", mDimName);
   iOptions.getValue("varName", mVarName);

   int dTime = getTimeDim(file);
   int dLon  = getLonDim(file);
   int dLat  = getLatDim(file);
   int dCoeff = getCoefficientDim(file);
   int nTime  = getDimSize(file, dTime);
   int nLat   = getDimSize(file, dLat);
   int nLon   = getDimSize(file, dLon);
   int nCoeff = getDimSize(file, dCoeff);


   // TODO: Deal with case where lat and lon are one-dimensional variables
   // Get latitudes
   int vLat = Util::MV;
   if(hasVar(file, "lat"))
      vLat = getVar(file, "lat");
   else if(hasVar(file, "latitude"))
      vLat = getVar(file, "latitude");
   else
      Util::error("Could not determine latitude variable");
   vec2 lats = getGridValues(file, vLat);

   // Get longitudes
   int vLon = Util::MV;
   if(hasVar(file, "lon"))
      vLon = getVar(file, "lon");
   else if(hasVar(file, "longitude"))
      vLon = getVar(file, "longitude");
   else
      Util::error("Could not determine longitude variable");
   vec2 lons = getGridValues(file, vLon);

   // Get elevations
   vec2 elevs;
   if(hasVar(file, "altitude")) {
      int vElev = getVar(file, "altitude");
      elevs = getGridValues(file, vElev);
   }
   else {
      elevs.resize(nLat);
      for(int i = 0; i < nLat; i++) {
         elevs[i].resize(nLon);
         for(int j = 0; j < nLon; j++) {
            elevs[i][j]= Util::MV;
         }
      }
   }

   int vTime = getVar(file, "time");
   double* times = new double[nTime];
   status = nc_get_var_double(file, vTime, times);
   handleNetcdfError(status, "could not get times");

   int var = getVar(file, mVarName);
   int totalNumParameters = nLat*nLon*nTime*nCoeff;
   float* values = getNcFloats(file, var);

   // Intitialize parameters to empty and then fill in later
   std::vector<float> params(nCoeff, Util::MV);
   Parameters parameters(params);
   for(int t = 0; t < nTime; t++) {
      mTimes.push_back(t);
      int time = t;
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int locIndex = i * nLon + j;
            Location location(lats[i][j], lons[i][j], elevs[i][j]);
            mParameters[location][time] = parameters;
         }
      }
   }

   /*
   Read parameters from file and arrange them in mParameters. This is a bit tricky because we do not
   force a certain ordering of dimensions in the coefficients variable. In addition, the variable
   may or may not have a time dimension. The algorithm is to figure out what order the various
   dimensions are and then loop over the retrieved parameters placing them into the right position
   in mParameters.
   */

   // Which index in list of dimension is each dimension?
   int latDimIndex = Util::MV;
   int lonDimIndex = Util::MV;
   int timeDimIndex = Util::MV;
   int coeffDimIndex = Util::MV;

   int ndims;
   status = nc_inq_varndims(file, var, &ndims);
   handleNetcdfError(status, "could not get number of dimensions of coefficients variable");

   // Fetch the dimensions of the coefficients variable
   int dims[ndims];
   status = nc_inq_vardimid(file, var, dims);
   handleNetcdfError(status, "could not get dimensions for coefficients variable");

   // Get dimensions sizes
   std::vector<int> sizes;
   for(int d = 0; d < ndims; d++) {
      size_t size;
      int status = nc_inq_dimlen(file, dims[d], &size);
      handleNetcdfError(status, "could not get size of a dimension for coefficients variable");
      sizes.push_back(size);
      if(dims[d] == dLat)
         latDimIndex = d;
      else if(dims[d] == dLon)
         lonDimIndex = d;
      else if(dims[d] == dTime)
         timeDimIndex = d;
      else if(dims[d] == dCoeff)
         coeffDimIndex = d;
   }

   // Check that coefficients has all required dimensions (does not need time)
   if(!Util::isValid(latDimIndex))
      Util::error("Coefficients in " + getFilename() + " is missing latitude dimension");
   if(!Util::isValid(lonDimIndex))
      Util::error("Coefficients in " + getFilename() + " is missing longitude dimension");
   if(!Util::isValid(coeffDimIndex))
      Util::error("Coefficients in " + getFilename() + " is missing coefficient dimension");

   // Loop over the parameters placing them into the right position
   int index = 0;
   while(index < totalNumParameters) {
      float currParameter = values[index];
      // Translate the linear index into a vector index
      std::vector<int> indices = getIndices(index, sizes);
      assert(indices.size() == ndims);

      // Determine the location, time, and parameter corresponding to 'index'
      int i = indices[latDimIndex];
      int j = indices[lonDimIndex];
      int coeffIndex = indices[coeffDimIndex];
      int timeIndex = 0;
      if(Util::isValid(timeDimIndex))
         timeIndex = indices[timeDimIndex];

      Location location(lats[i][j], lons[i][j], elevs[i][j]);

      // Assign parameter
      mParameters[location][timeIndex][coeffIndex] = currParameter;
      index++;
   }

   delete[] values;
   delete[] times;
}

int ParameterFileNetcdf::getDim(int iFile, std::string iDim) const {
   int dim;
   int status = nc_inq_dimid(iFile, iDim.c_str(), &dim);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return dim;
}

int ParameterFileNetcdf::getVar(int iFile, std::string iVar) const {
   int var;
   int status = nc_inq_varid(iFile, iVar.c_str(), &var);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return var;
}

bool ParameterFileNetcdf::hasDim(int iFile, std::string iDim) {
   int dim;
   int status = nc_inq_dimid(iFile, iDim.c_str(), &dim);
   return status == NC_NOERR;
}

bool ParameterFileNetcdf::hasVar(int iFile, std::string iVar) {
   int var;
   int status = nc_inq_varid(iFile, iVar.c_str(), &var);
   return status == NC_NOERR;
}

int ParameterFileNetcdf::getDimSize(int iFile, int iDim) const {
   size_t len;
   int status = nc_inq_dimlen(iFile, iDim, &len);
   return len;
}

std::vector<int> ParameterFileNetcdf::getDims(int iFile, int iVar) const {
   int ndims;
   int status = nc_inq_varndims(iFile, iVar, &ndims);
   handleNetcdfError(status, "could not get number of dimensions of coefficients variable");

   // Fetch the dimensions of the coefficients variable
   int dims[ndims];
   status = nc_inq_vardimid(iFile, iVar, dims);
   handleNetcdfError(status, "could not get dimensions for coefficients variable");
   std::vector<int> dimsVec(dims, dims+ndims);
   return dimsVec;
}

std::vector<int> ParameterFileNetcdf::getTimes() const {
   return mTimes;
}

bool ParameterFileNetcdf::isValid(std::string iFilename) {
   return true;
}

std::string ParameterFileNetcdf::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-p netcdf", "Parameters stored in a Netcdf file. File must have contain: dimensions time, lat (or latitude or y), lon (or longitude or x), coefficient; variables with dims: time[time], latitude[lat,lon], longitude[lat,lon], coefficients[*]. the coefficient variable must have lat, lon, coefficient dimensions, and optimally time. These can be in any order. If there is no time dimension, then the parameters are used for all times. The number of parameters in a set must be constant and equals the size of the 'coefficient' dimension.") << std::endl;
   ss << Util::formatDescription("   dimName=coefficient", "What is the name of the dimension representing different coefficients?") << std::endl;
   ss << Util::formatDescription("   varName=coefficients", "What is the name of the variable containing the coefficients?") << std::endl;
   ss << Util::formatDescription("   file=required", "Filename of file.") << std::endl;
   return ss.str();
}
void ParameterFileNetcdf::write() const {
   std::vector<Location> locations = getLocations();
   std::vector<int> times = getTimes();
}

void ParameterFileNetcdf::handleNetcdfError(int status, std::string message) const {
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

int ParameterFileNetcdf::getLatDim(int iFile) const {
   int dLat;
   if(hasDim(iFile, "lat"))
      dLat = getDim(iFile, "lat");
   else if(hasDim(iFile, "latitude"))
      dLat = getDim(iFile, "latitude");
   else if(hasDim(iFile, "y"))
      dLat = getDim(iFile, "y");
   else
      Util::error("Could not determine latitude dimension in " + getFilename());
   return dLat;
}

int ParameterFileNetcdf::getLonDim(int iFile) const {
   int dLon;
   if(hasDim(iFile, "lon"))
      dLon = getDim(iFile, "lon");
   else if(hasDim(iFile, "longitude"))
      dLon = getDim(iFile, "longitude");
   else if(hasDim(iFile, "x"))
      dLon = getDim(iFile, "x");
   else
      Util::error("Could not determine longitude dimension in " + getFilename());
   return dLon;
}

int ParameterFileNetcdf::getTimeDim(int iFile) const {
   int dTime;
   if(hasDim(iFile, "time"))
      dTime = getDim(iFile, "time");
   else
      Util::error("Could not determine time dimension in " + getFilename());
   return dTime;
}

int ParameterFileNetcdf::getCoefficientDim(int iFile) const {
   int dCoeff;
   if(hasDim(iFile, mDimName))
      dCoeff = getDim(iFile,mDimName);
   else
      Util::error("Could not find the coefficient dimension in " + getFilename());
   return dCoeff;
}

std::vector<int> ParameterFileNetcdf::getIndices(int i, const std::vector<int>& iCount) const {
   // The last index changes fastest, the first slowest
   int numDims = iCount.size();
   std::vector<int> indices;
   indices.resize(numDims);
   int sizeSoFar = 1;
   for(int k = numDims-1; k >= 0; k--) {
      int index = floor(i / sizeSoFar);
      index     = index % iCount[k];
      indices[k] = index;
      sizeSoFar *= iCount[k];
   }
   return indices;
}

float* ParameterFileNetcdf::getNcFloats(int iFile, int iVar) {
   int size = NetcdfUtil::getTotalSize(iFile, iVar);
   float* values = new float[size];
   int status = nc_get_var_float(iFile, iVar, values);
   char name[512];
   status = nc_inq_varname(iFile, iVar, name);
   handleNetcdfError(status, "could not get variable name");

   std::stringstream ss;
   ss <<  "could not retrieve values from variable '" << name << "'";
   handleNetcdfError(status, ss.str());

   // Convert missing values
   float MV = NetcdfUtil::getMissingValue(iFile, iVar);
   for(int i = 0; i < size; i++) {
      if(values[i] == MV)
         values[i] = Util::MV;
   }
   return values;
}

vec2 ParameterFileNetcdf::getGridValues(int iFile, int iVar) const {
   // Initialize values
   vec2 grid;
   int dLat = getLatDim(iFile);
   int dLon = getLonDim(iFile);
   int nLat = getDimSize(iFile, dLat);
   int nLon = getDimSize(iFile, dLon);
   grid.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      grid[i].resize(nLon, Util::MV);
   }

   // We have a lat/lon grid, where lat/lons are only provided along the pertinent dimension
   // Values are assumed to be constant across the other dimension.
   std::vector<int> dims = getDims(iFile, iVar);
   int numDims = dims.size();
   if(numDims == 1) {
      int dim = dims[0];
      long size = getDimSize(iFile, dim);
      float* values = new float[size];
      nc_get_var_float(iFile, iVar, values);
      // Latitude variable
      if(dim == getLatDim(iFile)) {
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               grid[i][j] = values[i];
            }
         }
      }
      // Longitude variable
      else if(dim == getLonDim(iFile)) {
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
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
      long count[numDims];
      int size = 1;
      int indexLat = Util::MV;
      int indexLon = Util::MV;
      std::vector<int> dims = getDims(iFile, iVar);
      for(int i = 0; i < numDims; i++) {
         if(dims[i] == getLatDim(iFile)) {
            count[i] = nLat;
            size *= count[i];
            indexLat = i;
         }
         else if(dims[i] == getLonDim(iFile)) {
            count[i] = nLon;
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
      nc_get_var_float(iFile, iVar, values);
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            // Latitude dimension is ordered first
            if(indexLat < indexLon) {
               grid[i][j] = values[i*nLon + j];
            }
            // Longitude dimension is ordered first
            else {
               grid[i][j] = values[j*nLat + i];
            }
         }
      }
      delete[] values;
   }
   return grid;
}
