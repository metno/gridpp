#include "MetnoKalman.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFileMetnoKalman::ParameterFileMetnoKalman(const Options& iOptions, bool iIsNew) : ParameterFile(iOptions, iIsNew),
      mLocalMV(-99999) {
   if(!isValid(getFilename()))
      Util::error(getFilename() + " is not a valid parameter file");

   int coeffFreq = 3; // How often are the coefficients for? In time steps.

   std::ifstream ifs(getFilename().c_str(), std::ifstream::in);

   // Empty file
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
      assert(status);

      for(int t = 0; t < (numTimes-1)*coeffFreq+1; t++) {
         mTimes.push_back(t);
      }
   }
   while(ifs.good()) {
      ifs.getline(line, 10000, '\n');
      if(ifs.good() && line[0] != '#') {
         std::stringstream ss(line);
         int stationId;
         bool status = ss >> stationId;
         assert(status);
         translate(stationId);

         float lat;
         status = ss >> lat;
         assert(status);
         translate(lat);

         float lon;
         status = ss >> lon;
         assert(status);
         translate(lon);

         float elev;
         status = ss >> elev;
         assert(status);
         translate(elev);

         float modelElev;
         status = ss >> modelElev;
         assert(status);
         translate(modelElev);

         float lastObs;
         status = ss >> lastObs;
         assert(status);
         translate(lastObs);

         // Loop over each value
         std::vector<float> values;
         while(ss.good()) {
            float value;
            bool status  = ss >> value;
            assert(status);
            translate(value);
            values.push_back(value);
         }
         assert(values.size() == numTimes);

         Location location(lat,lon,elev);
         // Values are for every 3 (coeffFreq) hours. Interpolate between so we get values every hour
         for(int i = 0; i < values.size(); i++) {
            setParameters(Parameters(values[i]), i*coeffFreq, location);
            // Fill in the gaps
            if(i < values.size() - 1) {
               for(int k = 1; k < coeffFreq; k++) {
                  float a = 1 - ((float) k)/coeffFreq;
                  float b = 1 - a;
                  if(Util::isValid(values[i]) && Util::isValid(values[i+1])) {
                     setParameters(Parameters(a*values[i]+b*values[i+1]), i*coeffFreq+k, location);
                  }
                  else {
                     setParameters(Parameters(Util::MV), i*coeffFreq+k, location);
                  }
               }
            }
         }
      }
   }
   ifs.close();

   recomputeTree();
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
   std::stringstream ss2(line);
   while(ss2 >> value) {
      values.push_back(value);
   }
   if(values.size() != 2)
      return false;
   int numTimes = values[1];
   int numCols = 6 + numTimes;
   while(true) {
      values.clear();
      ifs.getline(line, 10000, '\n');
      if(ifs.good() && line[0] != '#') {
         std::stringstream ss(line);
         float value;
         while(true) {
            if(!ss.good())
               break;
            bool status = ss >> value;
            values.push_back(value);
            if(!status)
               return false;
         }
         if(values.size() != numCols) {
            return false;
         }
      }
      else {
         break;
      }
   }
   return true;
}
bool ParameterFileMetnoKalman::isReadable() const {
   return ParameterFileMetnoKalman::isValid(getFilename());
}

std::string ParameterFileMetnoKalman::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-p metnoKalman", "Metno's internal format for kalman filter corrections.") << std::endl;
   ss << Util::formatDescription("   file=required", "Filename of file.") << std::endl;
   return ss.str();
}
