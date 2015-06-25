#include "ParameterFile.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFileMetnoKalman::ParameterFileMetnoKalman(std::string iFilename) : 
      mFilename(iFilename),
      mLocalMV(-99999) {
   std::ifstream ifs(mFilename.c_str(), std::ifstream::in);
   int coeffFreq = 3; // How often are the coefficients for? In time steps.

   if(!ifs.good()) {
      return;
   }
   std::vector<std::pair<std::pair<int, Location>, Parameters> > allParameters; // Time, Location, Parameters
   char line[10000];
   ifs.getline(line, 10000, '\n'); // Header line
   ifs.getline(line, 10000, '\n'); // Header line
   int temp     = Util::MV;
   int numTimes = Util::MV;
   if(ifs.good() && line[0] != '#') {
      std::stringstream ss(line);
      // Note: The first value read is not the number of locations. The true number is a bit fewer.
      bool status = ss >> temp >> numTimes;
      if(!status) {
         Util::error("Could not read second header line from file '" + mFilename + "'");
      }
   }
   while(ifs.good()) {
      ifs.getline(line, 10000, '\n');
      if(ifs.good() && line[0] != '#') {
         std::stringstream ss(line);
         int stationId;
         bool status = ss >> stationId;
         if(!status) {
            Util::error("Could not read station ID from file '" + mFilename + "'");
         }

         float lat;
         status = ss >> lat;
         if(!status) {
            Util::error("Could not read lat from file '" + mFilename + "'");
         }

         float lon;
         status = ss >> lon;
         if(!status) {
            Util::error("Could not read lon from file '" + mFilename + "'");
         }

         float elev;
         status = ss >> elev;
         if(!status) {
            Util::error("Could not read elev from file '" + mFilename + "'");
         }

         float modelElev;
         status = ss >> modelElev;
         if(!status) {
            Util::error("Could not read model elevation from file '" + mFilename + "'");
         }

         float lastObs;
         status = ss >> lastObs;
         if(!status) {
            Util::error("Could not read last observation from file '" + mFilename + "'");
         }
         if(lastObs == mLocalMV)
            lastObs = Util::MV;

         // Loop over each value
         std::vector<float> values;
         while(ss.good()) {
            float value;
            bool status  = ss >> value;
            if(!status) {
               Util::error("Could not read value from file '" + mFilename + "'");
            }
            if(value == mLocalMV)
               value = Util::MV;
            values.push_back(value);
         }
         if(values.size() != numTimes) {
            std::stringstream ss;
            ss << "ParameterFileMetnoKalman: Expecting " << numTimes << " columns of bias, found " << values.size();
            Util::error(ss.str());
         }
         Location location(lat,lon,elev);
         // Values are for every 3 (coeffFreq) hours. Interpolate between so we get values every hour
         std::vector<float> allValues(coeffFreq*(numTimes-1)+1, Util::MV);
         for(int i = 0; i < values.size(); i++) {
            allValues[i*coeffFreq]   = values[i];
            mTimes.push_back(i*coeffFreq);
            // Fill in the gaps
            if(i < values.size() - 1) {
               for(int k = 1; k < coeffFreq; k++) {
                  float a = 1 - ((float) k)/coeffFreq;
                  float b = 1 - a;
                  if(Util::isValid(values[i]) && Util::isValid(values[i+1])) {
                     allValues[i*coeffFreq+k] = (a*values[i]+b*values[i+1]);
                  }
                  mTimes.push_back(i*coeffFreq+k);
                  // allValues[i*3+1] = (values[i]*2+values[i+1]*1)/3;
                  // allValues[i*3+2] = (values[i]*1+values[i+1]*2)/3;
                  // mTimes.push_back(i*3+1);
                  // mTimes.push_back(i*3+2);
               }
            }
         }
         mParameters[location] = allValues;
      }
   }
   ifs.close();
}

Parameters ParameterFileMetnoKalman::getParameters(int iTime, const Location& iLocation) const {
   std::map<Location, std::vector<float> >::const_iterator it = mParameters.find(iLocation);
   if(it == mParameters.end())
      return Parameters();
   if(it->second.size() <= iTime) {
      std::stringstream ss;
      ss << "ParameterFileMetnoKalman: Retriving parameters for time " << iTime << ".";
      Util::error(ss.str());
   }
   float value = it->second[iTime];
   std::vector<float> param(1, value);
   return Parameters(param);
}

std::vector<Location> ParameterFileMetnoKalman::getLocations() const {
   std::vector<Location> locations;
   std::map<Location, std::vector<float> >::const_iterator it;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      locations.push_back(it->first);
   }
   return locations;
}

std::vector<int> ParameterFileMetnoKalman::getTimes() const {
   return mTimes;
}

std::string ParameterFileMetnoKalman::getFilename() const {
   return mFilename;
}

void ParameterFileMetnoKalman::setParameters(float iValue, int iTime, const Location& iLocation) {
   std::map<Location, std::vector<float> >::iterator it = mParameters.find(iLocation);
   if(iTime > mTimes[mTimes.size()-1]) {
      std::stringstream ss;
      ss << "ParameterFileMetnoKalman: Could not write parameters to time " << iTime;
      Util::error(ss.str());
   }
   if(it == mParameters.end()) {
      mParameters[iLocation].resize(mTimes.size(), Util::MV);
   }
   mParameters[iLocation][iTime] = iValue;
}
