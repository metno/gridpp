#include "Text.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>
#include <algorithm>

ParameterFileText::ParameterFileText(const Options& iOptions) : ParameterFile(iOptions),
      mIsSpatial(false) {
   std::ifstream ifs(getFilename().c_str(), std::ifstream::in);
   iOptions.getValue("spatial", mIsSpatial);
   if(!ifs.good()) {
      return;
   }
   mNumParameters = Util::MV;
   int counter = 0;
   std::set<int> times;
   while(ifs.good()) {
      char line[10000];
      ifs.getline(line, 10000, '\n');
      if(ifs.good() && line[0] != '#') {
         std::stringstream ss(line);
         // Loop over each value
         std::vector<float> values;
         int time;
         bool status = ss >> time;
         if(!status) {
            Util::error("Could not read time from file '" + mFilename + "'");
         }
         times.insert(time);

         Location location(Util::MV, Util::MV, Util::MV);
         if(mIsSpatial) {
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
            location = Location(lat, lon, elev);
         }

         while(ss.good()) {
            float value;
            bool status  = ss >> value;
            if(!status) {
               Util::error("Could not read value from file '" + mFilename + "'");
            }
            values.push_back(value);
         }
         if(mNumParameters == Util::MV)
            mNumParameters = values.size();
         else if(values.size() != mNumParameters) {
            std::stringstream ss;
            ss << "Parameter file '" + getFilename() + "' is corrupt, because it does not have the same"
               << " number of columns on each line" << std::endl;
            Util::error(ss.str());
         }
         Parameters parameters(values);
         mParameters[location][time] = parameters;
         counter++;
      }
   }
   ifs.close();
   mTimes = std::vector<int>(times.begin(), times.end());
   std::sort(mTimes.begin(), mTimes.end());

   // Ensure all locations have parameters for all times
   // Inserting empty parameter sets where there are none
   std::map<Location, std::map<int, Parameters> >::iterator it;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      std::map<int, Parameters>::const_iterator it2;
      for(int t = 0; t < mTimes.size(); t++) {
         if(it->second.find(t) == it->second.end()) {
            it->second[t] = Parameters();
         }
      }
   }
   std::stringstream ss;
   ss << "Reading " << mFilename << ". Found " << counter << " parameter sets.";
   Util::status(ss.str());
   if(!Util::isValid(mNumParameters))
      mNumParameters = 0;
}

std::vector<int> ParameterFileText::getTimes() const {
   return mTimes;
}

bool ParameterFileText::isFixedSize() const {
   return true;
}

bool ParameterFileText::isValid(std::string iFilename) {
   std::ifstream ifs(iFilename.c_str(), std::ifstream::in);
   if(!ifs.good()) {
      return false;
   }
   return true;
}

bool ParameterFileText::isReadable() const {
   std::ifstream ifs(getFilename().c_str(), std::ifstream::in);
   if(!ifs.good()) {
      return false;
   }
   return true;
}

int ParameterFileText::getNumParameters() const {
   return mNumParameters;
}

void ParameterFileText::write() const {
   write(mFilename);
}
void ParameterFileText::write(const std::string& iFilename) const {
   std::string filename = iFilename;
   std::ofstream ofs(filename.c_str(), std::ios_base::trunc);
   if(!ofs.good()) {
      Util::error("Cannot write parameters to " + filename);
   }

   std::map<Location, std::map<int, Parameters> >::const_iterator it;
   // Loop over times
   if(mIsSpatial)
      ofs << "# time lat lon elev parameters" << std::endl;
   else
      ofs << "# time parameters" << std::endl;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      std::map<int, Parameters>::const_iterator it2;
      // Loop over locations
      const Location& location = it->first;
      for(it2 = it->second.begin(); it2 != it->second.end(); it2++) {
         int time = it2->first;
         const Parameters& parameters = it2->second;
         if(parameters.size() != 0) {
            ofs << time;
            if(mIsSpatial) {
               ofs << " " << location.lat() << " " << location.lon() << " " << location.elev();
            }
            // Loop over parameter values
            for(int i = 0; i < parameters.size(); i++) {
               ofs << " " << parameters[i];
            }
            ofs << std::endl;
         }
      }
   }
   ofs.close();
}

bool ParameterFileText::isLocationDependent() const {
   return mIsSpatial;
}

std::string ParameterFileText::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-p text", "Simple ascii text file with space separated entries. Each line represents one lead time, and each column one parameter.") << std::endl;
   ss << Util::formatDescription("", "offset1 a b ... N") << std::endl;
   ss << Util::formatDescription("", "...") << std::endl;
   ss << Util::formatDescription("", "offsetM a b ... N") << std::endl;
   ss << Util::formatDescription("   file=required", "Filename of file.") << std::endl;
   ss << Util::formatDescription("   spatial=0", "Does this file contain spatial information? If yes, then the each line represents a leadtime and location (order of lines is not important):") << std::endl;
   ss << Util::formatDescription("", "offset1 lat1 lon1 elev1 a b ... N") << std::endl;
   ss << Util::formatDescription("", "                ...") << std::endl;
   ss << Util::formatDescription("", "offsetM lat1 lon1 elev1 a b ... N") << std::endl;
   ss << Util::formatDescription("", "                ...") << std::endl;
   ss << Util::formatDescription("", "offset1 latP lonP elevP a b ... N") << std::endl;
   ss << Util::formatDescription("", "                ...") << std::endl;
   ss << Util::formatDescription("", "offsetM latP lonP elevP a b ... N") << std::endl;
   return ss.str();
}
