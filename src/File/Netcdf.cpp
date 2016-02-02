#include "Netcdf.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileNetcdf::FileNetcdf(std::string iFilename, const Options& iOptions, bool iReadOnly) :
      File(iFilename, iOptions),
      mInDataMode(true) {
   int status = nc_open(getFilename().c_str(), iReadOnly ? NC_NOWRITE: NC_WRITE, &mFile);
   if(status != NC_NOERR) {
      Util::error("Could not open NetCDF file " + getFilename());
   }
}

FileNetcdf::~FileNetcdf() {
   nc_close(mFile);
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

int FileNetcdf::getDim(std::string iDim) const {
   int dim;
   int status = nc_inq_dimid(mFile, iDim.c_str(), &dim);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return dim;
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
   size_t len;
   int status = nc_inq_dimlen(mFile, dim, &len);
   return len;
}

int FileNetcdf::getDimSize(int iDim) const {
   size_t len;
   int status = nc_inq_dimlen(mFile, iDim, &len);
   return len;
}

int FileNetcdf::getNumDims(int iVar) const {
   int len;
   int status = nc_inq_varndims(mFile, iVar, &len);
   return len;
}

bool FileNetcdf::hasDim(std::string iDim) const {
   return hasDim(mFile, iDim);
}
bool FileNetcdf::hasVar(std::string iVar) const {
   return hasVar(mFile, iVar);
}

bool FileNetcdf::hasDim(int iFile, std::string iDim) {
   int dim;
   int status = nc_inq_dimid(iFile, iDim.c_str(), &dim);
   return status == NC_NOERR;
}
bool FileNetcdf::hasVar(int iFile, std::string iVar) {
   int var;
   int status = nc_inq_varid(iFile, iVar.c_str(), &var);
   return status == NC_NOERR;
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
      int dTime  = getDim("time");
      int id;
      int status = nc_def_var(mFile, "time", NC_DOUBLE, 1, &dTime, &id);
      handleNetcdfError(status, "creating time variable");
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
