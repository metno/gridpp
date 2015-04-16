#include "Netcdf.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileNetcdf::FileNetcdf(std::string iFilename, bool iReadOnly) :
      File(iFilename), 
      mFile(NcFile(getFilename().c_str(), iReadOnly ? NcFile::ReadOnly : NcFile::Write)) {
   if(!mFile.is_valid()) {
      Util::error("Error: Netcdf file " + getFilename() + " not valid");
   }
}

FileNetcdf::~FileNetcdf() {
   mFile.close();
}

bool FileNetcdf::hasVariableCore(Variable::Type iVariable) const {
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

bool FileNetcdf::hasDim(std::string iDim) const {
   return hasDim(mFile, iDim);
}
bool FileNetcdf::hasVar(std::string iVar) const {
   return hasVar(mFile, iVar);
}

bool FileNetcdf::hasDim(const NcFile& iFile, std::string iDim) {
   NcError q(NcError::silent_nonfatal); 
   NcDim* dim = iFile.get_dim(iDim.c_str());
   return dim != NULL;
}
bool FileNetcdf::hasVar(const NcFile& iFile, std::string iVar) {
   NcError q(NcError::silent_nonfatal); 
   NcVar* var = iFile.get_var(iVar.c_str());
   return var != NULL;
}

float FileNetcdf::getMissingValue(const NcVar* iVar) {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* fillValueAtt = iVar->get_att("_FillValue");
   if(fillValueAtt != NULL)
      return fillValueAtt->as_float(0);
   else
      return ncBad_float;
}
void FileNetcdf::setMissingValue(NcVar* iVar, float iValue) {
   // TODO: Mysterious errors can occur if the existing file was written
   // using an older HDF5 implementation, and the variable has 8 or more
   // attributes already. For more information, see "Corruption Problem
   // In HDF5 1.8.0 through HDF5 1.8.4" on
   // http://www.hdfgroup.org/HDF5/release/known_problems/index.html
   if(iValue != ncBad_float)
      iVar->add_att("_FillValue", iValue);
}

void FileNetcdf::setAttribute(NcVar* iVar, std::string iName, std::string iValue) {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* att = iVar->get_att(iName.c_str());
   if(att != NULL) {
      att->remove();
   }
   iVar->add_att(iName.c_str(), iValue.c_str());
}

void FileNetcdf::setGlobalAttribute(std::string iName, std::string iValue) {
   NcError q(NcError::silent_nonfatal); 
   NcAtt* att = mFile.get_att(iName.c_str());
   if(att != NULL) {
      att->remove();
   }
   mFile.add_att(iName.c_str(), iValue.c_str());
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
         times[i] = ncBad_double;
   }
   if(!hasVar("time")) {
      NcDim* dTime    = getDim("time");
      mFile.add_var("time", ncDouble, dTime);
   }
   NcVar* vTime = getVar("time");
   double timesArr[getNumTime()];
   for(int t = 0; t < getNumTime(); t++) {
      timesArr[t] = times[t];
   }
   vTime->put(timesArr, getNumTime());
   setAttribute(vTime, "long_name", "time");
   setAttribute(vTime, "standard_name", "time");
   setAttribute(vTime, "units", "seconds since 1970-01-01 00:00:00 +00:00");
}
void FileNetcdf::writeReferenceTime() {
   if(!hasVar("forecast_reference_time")) {
      mFile.add_var("forecast_reference_time", ncDouble);
   }
   NcVar* vTime = getVar("forecast_reference_time");
   double referenceTime = getReferenceTime();
   if(!Util::isValid(referenceTime))
      referenceTime = ncBad_double;
   vTime->put(&referenceTime, 1);
   setAttribute(vTime, "standard_name", "forecast_reference_time");
   setAttribute(vTime, "units", "seconds since 1970-01-01 00:00:00 +00:00");
}
void FileNetcdf::writeGlobalAttributes() {
   setGlobalAttribute("Conventions", "CF-1.0");
   std::stringstream ss;
   ss << Util::getCurrentDateString() << ": creation by gridpp";
   setGlobalAttribute("history", ss.str());
}
