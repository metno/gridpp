#include "NorcomQnh.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include "../Util.h"
#include <iomanip>

FileNorcomQnh::FileNorcomQnh(std::string iFilename, const Options& iOptions) :
      File(iFilename, iOptions) {
   mLats.resize(1);
   mLons.resize(1);
   mElevs.resize(1);
   mLandFractions.resize(1);
   if(!iOptions.getValues("lats", mLats[0])) {
      Util::error("Missing 'lats' option for '" + iFilename + "'");
   }
   if(!iOptions.getValues("lons", mLons[0])) {
      Util::error("Missing 'lons' option for '" + iFilename + "'");
   }
   if(!iOptions.getValues("elevs", mElevs[0])) {
      Util::error("Missing 'elevs' option for '" + iFilename + "'");
   }
   mLandFractions[0].resize(mElevs[0].size(), Util::MV);
   if(!iOptions.getValues("names", mNames)) {
      Util::error("Missing 'names' option for '" + iFilename + "'");
   }
   if(!iOptions.getValue("numTimes", mNTime)) {
      Util::error("Missing 'numTimes' option for '" + iFilename + "'");
   }
   if(mLats[0].size() != mLons[0].size() || mLats[0].size() != mElevs[0].size() || mLats[0].size() != mNames.size()) {
      Util::error("FileNorcomQnh: 'lats', 'lons', 'elevs', 'names' must be the same size");
   }
   for(int i = 0; i < mLats[0].size(); i++) {
      float lat = mLats[0][i];
      if(lat < -90 || lat > 90) {
         std::stringstream ss;
         ss << "Invalid latitude: " << lat;
         Util::error(ss.str());
      }
   }
   mNLat = 1;
   mNLon = mLats[0].size();
   mNEns = 1;

   std::vector<double> times;
   for(int i = 0; i < mNTime; i++)
      times.push_back(i);

   // Determine the times for this filetype.
   if(!iOptions.getValue("startTime", mStartTime)) {
      Util::error("Missing 'startTime' option for '" + iFilename + "'");
   }
   if(!iOptions.getValue("endTime", mEndTime)) {
      Util::error("Missing 'endTime' option for '" + iFilename + "'");
   }
   if(mStartTime > mEndTime) {
      Util::error("FileNorcomQnh: 'startTime' must be <= 'endTime'");
   }
   setTimes(times);
}

FileNorcomQnh::~FileNorcomQnh() {
}

FieldPtr FileNorcomQnh::getFieldCore(Variable::Type iVariable, int iTime) const {
   FieldPtr field = getEmptyField();
   return field;
}

void FileNorcomQnh::writeCore(std::vector<Variable::Type> iVariables) {
   std::ofstream ofs(getFilename().c_str());
   if(iVariables.size() == 0) {
      Util::warning("No variables to write");
      return;
   }
   Variable::Type variable = iVariables[0];
   if(iVariables.size() > 1) {
      std::stringstream ss;
      ss <<"Output NorcomQnh can only write one variables, several given. Will write variable ";
      ss << Variable::getTypeName(variable);
      Util::warning(ss.str());
   }

   // Find the length of the longest station name
   int maxNameSize = 0;
   for(int j = 0; j < mNames.size(); j++) {
      if(mNames[j].size() > maxNameSize)
         maxNameSize = mNames[j].size();
   }

   // Write header
   time_t currUnixTime = Util::getCurrentUnixTime();
   std::string currTimeStamp = getNorcomTimeStamp(currUnixTime);
   std::vector<double> times = getTimes();
   std::string startTime = getNorcomTimeStamp(times[mStartTime]);
   std::string endTime   = getNorcomTimeStamp(times[mEndTime]);
   ofs << "FBNO52 ENNC " << currTimeStamp << std::endl;
   ofs << "VALID " << startTime << " - " << endTime << " UTC." << std::endl;

   ofs.precision(0);
   // Write one line for each station
   for(int j = 0; j < mLats[0].size(); j++) {
      std::string locationName = mNames[j];
      ofs << "EST MIN QNH ";
      ofs << std::setfill(' ') << std::setw(maxNameSize) << std::left << locationName << ": ";

      // Find minimum
      int valuePa = Util::MV;
      for(int t = mStartTime; t <= mEndTime; t++) {
         FieldPtr field = getField(variable, t);
         int currValue = (int) (*field)(0,j,0);
         if(!Util::isValid(valuePa)|| (Util::isValid(currValue) && currValue < valuePa))
            valuePa = currValue;
      }
      if(!Util::isValid(valuePa)) {
         Util::error("Invalid value when writing QNH to NorcomQnh");
      }
      int valueHpa = valuePa / 100;
      ofs << std::setfill('0') << std::setw(4) << std::right << valueHpa << " HPA" << std::endl;
   }
   ofs.close();
}

std::string FileNorcomQnh::getNorcomTimeStamp(time_t iUnixTime) const {
   int date = Util::getDate(iUnixTime);
   int time = Util::getTime(iUnixTime);
   std::stringstream ss;
   int day = date % 100; 
   int HHMM = time / 100;
   ss << std::setfill('0') << std::right << std::setw(2) << day;
   ss << std::setfill('0') << std::right << std::setw(4) << HHMM;

   return ss.str();
}

std::string FileNorcomQnh::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=norcomQnh", "Output format for sending minimum QNH values to Norcom. Finds the minimum QNH values on the interval [startTime,endTime]. QNH must either exist in the input file, or created using a calibrator (such as -c qnh)") << std::endl;
   ss << Util::formatDescription("   lats=required", "Comma-separated list of latitudes (lat1,lat2,lat3,...). Values in degrees, north is positive.") << std::endl;
   ss << Util::formatDescription("   lons=required", "Longitudes (in degrees, east is positive)") << std::endl;
   ss << Util::formatDescription("   elevs=required", "Elevations (in meters)") << std::endl;
   ss << Util::formatDescription("   names=required", "Comma-separated list of station names.") << std::endl;
   ss << Util::formatDescription("   numTimes=undef", "Number of times. Set this equal to the number of times in the input file.") << std::endl;
   ss << Util::formatDescription("   startTime=undef", "First time index to find minimum over.") << std::endl;
   ss << Util::formatDescription("   endTime=undef", "Last time index to find minimum over.") << std::endl;
   return ss.str();
}
