#include "Arome.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileArome::FileArome(std::string iFilename, bool iReadOnly) : FileNetcdf(iFilename, iReadOnly) {
   // Set dimensions
   NcDim* dTime = getDim("time");
   NcDim* dLon  = getDim("x");
   NcDim* dLat  = getDim("y");
   mNTime = dTime->size();
   mNLat  = dLat->size();
   mNLon  = dLon->size();
   mNEns  = 1;

   mLats = getLatLonVariable("latitude");
   mLons = getLatLonVariable("longitude");
   if(hasVariableCore("surface_geopotential")) {
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
   else {
      mElevs = getLatLonVariable("altitude");
   }

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

FieldPtr FileArome::getFieldCore(Variable::Type iVariable, int iTime) const {
   std::string variableName = getVariableName(iVariable);
   return getFieldCore(variableName, iTime);
}
FieldPtr FileArome::getFieldCore(std::string iVariable, int iTime) const {
   // Not cached, retrieve data
   NcVar* var = getVar(iVariable);
   int nLat  = mNLat;
   int nLon  = mNLon;

   int numDims = var->num_dims();

   long* count;
   long totalCount = nLat*nLon;
   if(numDims == 4) {
      // Variable has a surface dimension
      count = new long[4];
      count[0] = 1;
      count[1] = 1;
      count[2] = nLat;
      count[3] = nLon;
      var->set_cur(iTime, 0, 0, 0);
   }
   else if(numDims == 3) {
      count = new long[3];
      count[0] = 1;
      count[1] = nLat;
      count[2] = nLon;
      var->set_cur(iTime, 0, 0);
   }
   else {
      std::stringstream ss;
      ss << "Cannot read variable '" << iVariable << "' from '" << getFilename() << "'";
      Util::error(ss.str());
   }
   float* values = new float[nLat*nLon];
   var->get(values, count);
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
   return field;
}

FileArome::~FileArome() {
}

void FileArome::writeCore(std::vector<Variable::Type> iVariables) {
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
         if(0) {
            NcDim* dTime    = getDim("time");
            NcDim* dSurface = getDim("height0");
            NcDim* dLon     = getDim("x");
            NcDim* dLat     = getDim("y");
            var = mFile.add_var(variable.c_str(), ncFloat, dTime, dSurface, dLat, dLon);
         }
         else {
            NcDim* dTime    = getDim("time");
            NcDim* dLon     = getDim("x");
            NcDim* dLat     = getDim("y");
            var = mFile.add_var(variable.c_str(), ncFloat, dTime, dLat, dLon);
         }
      }
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
            int numDims = var->num_dims();
            if(numDims == 4) {
               var->set_cur(t, 0, 0, 0);
               var->put(values, 1, 1, mNLat, mNLon);
            }
            else if(numDims == 3) {
               var->set_cur(t, 0, 0);
               var->put(values, 1, mNLat, mNLon);
            }
            else {
               std::stringstream ss;
               ss << "Cannot write variable '" << variable << "' from '" << getFilename() << "'";
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
   NcVar* var = getVar(iVar);
   float MV = getMissingValue(var);
   long count[2] = {getNumLat(), getNumLon()};
   float* values = new float[getNumLon()*getNumLat()];
   var->get(values, count);
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
   bool status = false;
   NcFile file = NcFile(iFilename.c_str(), NcFile::ReadOnly);
   if(file.is_valid()) {
      status = hasDim(file, "time") && hasDim(file, "x") && hasDim(file, "y") &&
               !hasDim(file, "ensemble_member") &&
               hasVar(file, "latitude") && hasVar(file, "longitude");
   }
   file.close();
   return status;
}

std::string FileArome::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=arome", "AROME file") << std::endl;
   return ss.str();
}
