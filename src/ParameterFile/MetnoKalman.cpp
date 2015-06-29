#include "MetnoKalman.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFileMetnoKalman::ParameterFileMetnoKalman(std::string iFilename) : ParameterFile(iFilename),
      mLocalMV(-99999) {
   std::ifstream ifs(iFilename.c_str(), std::ifstream::in);
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
         for(int i = 0; i < values.size(); i++) {
            mParameters[location][i*coeffFreq] = Parameters(values[i]);
            mTimes.push_back(i*coeffFreq);
            // Fill in the gaps
            if(i < values.size() - 1) {
               for(int k = 1; k < coeffFreq; k++) {
                  float a = 1 - ((float) k)/coeffFreq;
                  float b = 1 - a;
                  if(Util::isValid(values[i]) && Util::isValid(values[i+1])) {
                     mParameters[location][i*coeffFreq+k] = Parameters(a*values[i]+b*values[i+1]);
                  }
                  else {
                     mParameters[location][i*coeffFreq+k] = Parameters(Util::MV);
                  }
                  mTimes.push_back(i*coeffFreq+k);
               }
            }
         }
      }
   }
   ifs.close();
}

std::vector<int> ParameterFileMetnoKalman::getTimes() const {
   return mTimes;
}

bool ParameterFileMetnoKalman::isValid(std::string iFilename) {
   std::ifstream ifs(iFilename.c_str(), std::ifstream::in);
   if(!ifs.good()) {
      return false;
   }

   char line[10000];
   ifs.getline(line, 10000, '\n'); // Header line
   std::stringstream ss(line);
   int value;
   std::vector<int> values;
   while(ss >> value) {
      values.push_back(value);
   }
   if(values.size() != 4)
      return false;

   values.clear();
   ifs.getline(line, 10000, '\n'); // Header line
   while(ss >> value) {
      values.push_back(value);
   }
   if(values.size() != 2)
      return false;
   return true;
}
