#include "Netcdf.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFileNetcdf::ParameterFileNetcdf(const Options& iOptions) : ParameterFile(iOptions),
      mDimName("coeff"),
      mVarName("coefficient") {
   int file;
   int status = nc_open(getFilename().c_str(), NC_NOWRITE, &file);
   handleNetcdfError(status, "invalid parameter file");
   iOptions.getValue("dimName", mDimName);
   iOptions.getValue("varName", mVarName);

   int dTime = getDim(file, "time");
   int dLon;
   if(hasDim(file, "lon"))
      dLon = getDim(file, "lon");
   else if(hasDim(file, "longitude"))
      dLon = getDim(file, "longitude");
   else if(hasDim(file, "x"))
      dLon = getDim(file, "x");
   else {
      Util::error("Could not determine longitude dimension in " + getFilename());
   }
   int dLat;
   if(hasDim(file, "lat"))
      dLat = getDim(file, "lat");
   else if(hasDim(file, "latitude"))
      dLat = getDim(file, "latitude");
   else if(hasDim(file, "y"))
      dLat = getDim(file, "y");
   else {
      Util::error("Could not determine latitude dimension in " + getFilename());
   }
   int dCoeff = getDim(file, mDimName);
   int nTime  = getDimSize(file, dTime);
   int nLat   = getDimSize(file, dLat);
   int nLon   = getDimSize(file, dLon);
   int nCoeff = getDimSize(file, dCoeff);

   long count2[2] = {nLat, nLon};
   // Get latitudes
   int vLat = Util::MV;
   if(hasVar(file, "lat"))
      vLat = getVar(file, "lat");
   else if(hasVar(file, "latitude"))
      vLat = getVar(file, "latitude");
   else
      Util::error("Could not determine latitude variable");
   float* lats = new float[nLat*nLon];
   status = nc_get_var_float(file, vLat, lats);
   handleNetcdfError(status, "could not get latitudes");

   // Get longitudes
   int vLon = Util::MV;
   if(hasVar(file, "lon"))
      vLon = getVar(file, "lon");
   else if(hasVar(file, "longitude"))
      vLon = getVar(file, "longitude");
   else
      Util::error("Could not determine longitude variable");
   float* lons = new float[nLat*nLon];
   status = nc_get_var_float(file, vLon, lons);
   handleNetcdfError(status, "could not get longitudes");

   // Get elevations
   float* elevs = new float[nLat*nLon];
   if(hasVar(file, "altitude")) {
      int vElev = getVar(file, "altitude");
      status = nc_get_var_float(file, vElev, elevs);
      handleNetcdfError(status, "could not get altitudes");
   }
   else {
      for(int i = 0; i < nLat*nLon; i++) {
         elevs[i] = Util::MV;
      }
   }

   int vTime = getVar(file, "time");
   double* times = new double[nTime];
   status = nc_get_var_double(file, vTime, times);
   handleNetcdfError(status, "could not get times");

   int var = getVar(file, mVarName);
   float* values = new float[nLat*nLon*nTime*nCoeff];
   status = nc_get_var_float(file, var, values);
   handleNetcdfError(status, "could not get parameter values");
   int index = 0;
   for(int t = 0; t < nTime; t++) {
      mTimes.push_back(t);
      int time = t;
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int locIndex = i * nLon + j;
            Location location(lats[locIndex], lons[locIndex], elevs[locIndex]);
            std::vector<float> params(nCoeff, Util::MV);
            for(int p = 0; p < nCoeff; p++) {
               params[p] = values[index];
               index++;
            }
            Parameters parameters(params);
            mParameters[location][time] = parameters;
         }
      }
   }
   delete[] lats;
   delete[] lons;
   delete[] elevs;
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


std::vector<int> ParameterFileNetcdf::getTimes() const {
   return mTimes;
}

bool ParameterFileNetcdf::isValid(std::string iFilename) {
   return true;
}

std::string ParameterFileNetcdf::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-p netcdf", "Parameters stored in a Netcdf file. File must have contain: dimensions time, lat (or y), lon (or x), coeff; variables with dims: time[time], latitude[lat,lon], longitude[lat,lon], coefficients[time,lat,lon,coeff]. The number of parameters in a set must be constant and equals the size of the 'coeff' dimension.") << std::endl;
   ss << Util::formatDescription("   dimName=coefficient", "What is the name of the dimension representing different coefficients?") << std::endl;
   ss << Util::formatDescription("   varName=coeff", "What is the name of the variable containing the coefficients?") << std::endl;
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
