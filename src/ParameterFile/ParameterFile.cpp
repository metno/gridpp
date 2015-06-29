#include "ParameterFile.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFile::ParameterFile(std::string iFilename) : mFilename(iFilename) {

}

ParameterFile* ParameterFile::getScheme(std::string iFilename, const Options& iOptions) {
   std::string name;
   if(ParameterFileMetnoKalman::isValid(iFilename)) {
      name = "metnoKalman";
   }
   else if(ParameterFileText::isValid(iFilename)) {
      name = "text";
   }
   else {
      Util::error("Could not determine type for parameter file '" + iFilename + "'");
   }
   return getScheme(name, iFilename, iOptions);
}

ParameterFile* ParameterFile::getScheme(std::string iName, std::string iFilename, const Options& iOptions) {
   ParameterFile* p;
   if(iName == "metnoKalman") {
      p = new ParameterFileMetnoKalman(iFilename);
   }
   else if(iName == "text") {
      p = new ParameterFileText(iFilename);
   }
   else if(iName == "textSpatial") {
      p = new ParameterFileText(iFilename, true);
   }
   else {
      Util::error("Could not determine type for parameter file '" + iFilename + "'");
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
   // Find the right location to use
   std::map<Location, std::map<int, Parameters> >::const_iterator it;
   if(mParameters.size() == 1) {
      // One set of parameters for all locations
       it = mParameters.begin();
   }
   else {
      it = mParameters.find(iLocation);
      // Nearest neighbour
      // TODO
      if(it == mParameters.end()) {
      }
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

std::string ParameterFile::getDescription() {
   std::stringstream ss;
   ss << "What file type is the parameters in? One of 'text', 'textSpatial, and 'metnoKalman'.";
   return ss.str();
}
