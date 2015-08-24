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
   NcFile file(getFilename().c_str(), NcFile::ReadOnly);
   if(!file.is_valid()) {
      Util::error("Parameter file " + getFilename() + " not valid");
   }
   iOptions.getValue("dimName", mDimName);
   iOptions.getValue("varName", mVarName);

   NcDim* dTime = getDim(file, "time");
   NcDim* dLon;
   if(hasDim(file, "lon"))
      dLon = getDim(file, "lon");
   else if(hasDim(file, "x"))
      dLon = getDim(file, "x");
   else {
      Util::error("Could not determine longitude dimension in " + getFilename());
   }
   NcDim* dLat;
   if(hasDim(file, "lat"))
      dLat = getDim(file, "lat");
   else if(hasDim(file, "y"))
      dLat = getDim(file, "y");
   else {
      Util::error("Could not determine latitude dimension in " + getFilename());
   }
   NcDim* dCoeff = getDim(file, mDimName);
   int nTime  = dTime->size();
   int nLat   = dLat->size();
   int nLon   = dLon->size();
   int nCoeff = dCoeff->size();

   // get lats lons
   long count2[2] = {nLat, nLon};
   NcVar* vLat = getVar(file, "latitude");
   NcVar* vLon = getVar(file, "longitude");
   NcVar* vElev = getVar(file, "altitude");
   float* lats = new float[nLat*nLon];
   float* lons = new float[nLat*nLon];
   float* elevs = new float[nLat*nLon];
   vLat->get(lats, count2);
   vLon->get(lons, count2);
   if(hasVar(file, "altitude")) {
      vElev->get(elevs, count2);
   }
   else {
      for(int i = 0; i < nLat*nLon; i++) {
         elevs[i] = Util::MV;
      }
   }

   NcVar* vTime = getVar(file, "time");
   double* times = new double[nTime];
   vTime->get(times , nTime);

   NcVar* var = getVar(file, mVarName);
   long count[4] = {nTime, nLat, nLon, nCoeff};
   float* values = new float[nLat*nLon*nTime*nCoeff];
   var->get(values, count);
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

NcDim* ParameterFileNetcdf::getDim(const NcFile& iFile, std::string iDim) const {
   NcError q(NcError::silent_nonfatal); 
   NcDim* dim = iFile.get_dim(iDim.c_str());
   if(dim == NULL) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return dim;
}
NcVar* ParameterFileNetcdf::getVar(const NcFile& iFile, std::string iVar) const {
   NcError q(NcError::silent_nonfatal); 
   NcVar* var = iFile.get_var(iVar.c_str());
   if(var == NULL) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return var;
}

std::vector<int> ParameterFileNetcdf::getTimes() const {
   return mTimes;
}

bool ParameterFileNetcdf::isValid(std::string iFilename) {
   return true;
}
bool ParameterFileNetcdf::hasDim(const NcFile& iFile, std::string iDim) const {
   NcError q(NcError::silent_nonfatal); 
   NcDim* dim = iFile.get_dim(iDim.c_str());
   return dim != NULL;
}
bool ParameterFileNetcdf::hasVar(const NcFile& iFile, std::string iVar) const {
   NcError q(NcError::silent_nonfatal); 

   NcVar* var = iFile.get_var(iVar.c_str());
   return var != NULL;
}

std::string ParameterFileNetcdf::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-p netcdf", "Parameters stored in a Netcdf file. File must have contain: dimensions time, lat (or y), lon (or x), coeff; variables with dims: time[time], latitude[lat,lon], longitude[lat,lon], coefficients[time,lat,lon,coeff]. The number of parameters in a set must be constant and equals the size of the 'coeff' dimension.") << std::endl;
   ss << Util::formatDescription("   dimName=coefficient", "What is the name of the dimension representing different coefficients?") << std::endl;
   ss << Util::formatDescription("   varName=coeff", "What is the name of the variable containing the coefficients?") << std::endl;
   ss << Util::formatDescription("   file=required", "Filename of file.") << std::endl;
   return ss.str();
}
