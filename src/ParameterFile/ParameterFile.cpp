#include "ParameterFile.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFile::ParameterFile(const Options& iOptions) : mFilename("") {
   iOptions.getValue("file", mFilename);
}

ParameterFile* ParameterFile::getScheme(std::string iName, const Options& iOptions) {
   ParameterFile* p;
   if(iName == "metnoKalman") {
      p = new ParameterFileMetnoKalman(iOptions);
   }
   else if(iName == "text") {
      p = new ParameterFileText(iOptions);
   }
   else if(iName == "textSpatial") {
      p = new ParameterFileText(iOptions, true);
   }
   else {
      Util::error("Parameter file type '" + iName + "' not recognized");
   }
   return p;
}

Parameters ParameterFile::getParameters(int iTime) const {
   // Find the right location to use
   std::map<Location, std::map<int, Parameters> >::const_iterator it;
   if(mParameters.size() == 1) {
      // One set of parameters for all locations
       it = mParameters.begin();
   }

   // Find the right time to use
   std::map<int,Parameters> timeParameters = it->second;
   if(timeParameters.size() == 1) {
      // One set of parameters for all times
      return timeParameters.begin()->second;
   }
   else if(it->second.find(iTime) != it->second.end()) {
      return timeParameters[iTime];
   }
   else {
      std::stringstream ss;
      ss << "Parameter file '" << mFilename << "' does not have values for time " << iTime << ".";
      Util::error(ss.str());
   }
}
Parameters ParameterFile::getParameters(int iTime, const Location& iLocation) const {
   if(mParameters.size() == 0)
      return Parameters();
   // Find the right location to use
   Location loc(Util::MV,Util::MV,Util::MV);
   if(mParameters.size() == 1) {
      // One set of parameters for all locations
      std::map<Location, std::map<int, Parameters> >::const_iterator it = mParameters.begin();
      loc = it->first;
   }
   else {
      std::map<Location, std::map<int, Parameters> >::const_iterator it;
      it = mParameters.find(iLocation);
      // Nearest neighbour
      // TODO
      if(it == mParameters.end()) {
         std::map<Location, std::map<int, Parameters> >::const_iterator it;
         float minDist = Util::MV;
         for(it = mParameters.begin(); it != mParameters.end(); it++) {
            const Location& currLoc = it->first;
            float dist = iLocation.getDistance(currLoc);
            if(!Util::isValid(minDist) || (Util::isValid(dist) && dist < minDist)) {
               minDist = dist;
               loc = it->first;
               break;
            }
         }
         if(!Util::isValid(minDist)) {
            // TODO: No nearby location found
            Util::error("No parameter locations found");
         }
         else {
         }
      }
      else {
         loc = iLocation;
      }
   }

   // Find the right time to use
   std::map<int,Parameters> timeParameters = mParameters.find(loc)->second;//it->second;
   if(timeParameters.size() == 1) {
      // One set of parameters for all times
      return timeParameters.begin()->second;
   }
   else if(timeParameters.find(iTime) != timeParameters.end()) {
      return timeParameters[iTime];
   }
   else {
      std::stringstream ss;
      ss << "Parameter file '" << mFilename << "' does not have values for time " << iTime << ".";
      Util::error(ss.str());
   }
}

void ParameterFile::setParameters(Parameters iParameters, int iTime, const Location& iLocation) {
   mParameters[iLocation][iTime] = iParameters;
}
void ParameterFile::setParameters(Parameters iParameters, int iTime) {
   std::map<Location, std::map<int,Parameters> >::iterator it = mParameters.begin();
   it->second[iTime] = iParameters;
}

std::string ParameterFile::getFilename() const {
   return mFilename;
}

std::vector<Location> ParameterFile::getLocations() const {
   std::vector<Location> locations;
   std::map<Location, std::map<int, Parameters> >::const_iterator it;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      locations.push_back(it->first);
   }
   return locations;
}

bool ParameterFile::isLocationDependent() const {
   return mParameters.size() > 1;
}
// TODO
bool ParameterFile::isTimeDependent() const {
   return true;
}

int ParameterFile::getNumParameters() const {
   int size = Util::MV;

   std::vector<Location> locations;
   std::map<Location, std::map<int, Parameters> >::const_iterator it;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      std::map<int, Parameters>::const_iterator it2;
      for(it2 = it->second.begin(); it2 != it->second.end(); it2++) {
         int currSize = it2->second.size();
         if(Util::isValid(size) && currSize != currSize)
            return Util::MV;
         size = currSize;
      }
   }
   return size;
}

std::string ParameterFile::getDescription(bool iSpatialOnly) {
   std::stringstream ss;
   if(iSpatialOnly)
      ss << "What file type is the parameters in? One of 'text' and 'metnoKalman'.";
   else
      ss << "What file type is the parameters in? One of 'text' and 'metnoKalman'.";
   return ss.str();
}

void ParameterFile::setFilename(std::string iFilename) {
   mFilename = iFilename;
}
