#include "Point.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include "../Util.h"

FilePoint::FilePoint(std::string iFilename, const Options& iOptions) :
      File(iFilename) {
   float lat,lon,elev;
   if(!iOptions.getValue("lat", lat)) {
      Util::error("Missing 'lat' option for '" + iFilename + "'");
   }
   if(!iOptions.getValue("lon", lon)) {
      Util::error("Missing 'lon' option for '" + iFilename + "'");
   }
   if(!iOptions.getValue("elev", elev)) {
      Util::error("Missing 'elev' option for '" + iFilename + "'");
   }
   if(lat < -90 || lat > 90) {
      std::stringstream ss;
      ss << "Invalid latitude: " << lat;
      Util::error(ss.str());
   }
   std::vector<float> lat0(1, lat);
   std::vector<float> lon0(1, lon);
   std::vector<float> elev0(1, elev);
   mLats.push_back(lat0);
   mLons.push_back(lon0);
   mElevs.push_back(elev0);
   mNLat = 1;
   mNLon = 1;
   mNEns = 1;
   std::vector<double> times;

   // Determine the times for this filetype.
   if(iOptions.getValue("time", mNTime)) {
      // Empty file, probably used as output only
      for(int i = 0; i < mNTime; i++)
         times.push_back(Util::MV);
   }
   else {
      // Existing file
      std::ifstream ifs(getFilename().c_str());
      if(ifs.good()) {
         while(ifs.good()) {
            double time, value;
            ifs >> time >> value;
            if(ifs.good())
               times.push_back(time);
         }
         mNTime = times.size();
      }
      else {
         Util::error("Missing 'time' option for empty file '" + iFilename + "'");
      }
   }
   setTimes(times);
}

FilePoint::~FilePoint() {
}

FieldPtr FilePoint::getFieldCore(Variable::Type iVariable, int iTime) const {
   std::ifstream ifs(getFilename().c_str());
   FieldPtr field = getEmptyField();
   for(int i = 0; i < getNumTime(); i++) {
      float time, value;
      ifs >> time >> value;
      if(i == iTime)
         (*field)(0,0,0) = value;
   }
   ifs.close();
   return field;
}

void FilePoint::writeCore(std::vector<Variable::Type> iVariables) {
   std::ofstream ofs(getFilename().c_str());
   // ofs << mLats[0][0] << " " << mLons[0][0] << " " << mElevs[0][0];
   if(iVariables.size() == 0) {
      Util::warning("No variables to write");
      return;
   }
   for(int i = 0; i < getNumTime(); i++) {
      ofs.precision(0);
      FieldPtr field = getField(iVariables[0], i);
      if(field != NULL) {
         ofs << (long) getTimes()[i];
         ofs.precision(2);
         ofs << std::fixed << " " << (*field)(0,0,0);
         ofs << std::endl;
      }
   }
   ofs.close();
}

std::string FilePoint::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=point", "Point file for one location. Each row is a time when the first column is the UNIX time and the second the forecast value.") << std::endl;
   ss << Util::formatDescription("   lat=required", "Latitude (in degrees, north is positive)") << std::endl;
   ss << Util::formatDescription("   lon=required", "Longitude (in degrees, east is positive)") << std::endl;
   ss << Util::formatDescription("   elev=required", "Elevation (in meters)") << std::endl;
   ss << Util::formatDescription("   time=undef", "Number of times. Required if the file does not exist.") << std::endl;
   return ss.str();
}
