#include "ParameterFile.h"
#include <fstream>
#include <sstream>
#include "Util.h"
#include <assert.h>

ParameterFile::ParameterFile(std::string iFilename) : 
      mFilename(iFilename) {
   std::ifstream ifs(mFilename.c_str(), std::ifstream::in);
   if(!ifs.good()) {
      Util::error("Parameter file '" + iFilename + "' does not exist");
   }
   int currSize = Util::MV;
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
         while(ss.good()) {
            float value;
            bool status  = ss >> value;
            if(!status) {
               Util::error("Could not read value from file '" + mFilename + "'");
            }
            values.push_back(value);
         }
         if(Util::isValid(currSize) && values.size() != currSize) {
            std::stringstream ss;
            ss << "Parameter file '" + iFilename + "' is corrupt, because it does not have the same"
               << " number of oclumns on each line" << std::endl;
            Util::error(ss.str());
         }
         Parameters parameters(values);
         mParameters[time] = parameters;
      }
   }
   ifs.close();
   std::stringstream ss;
   ss << "Reading " << mFilename << ". Found " << getSize() << " parameter sets.";
   Util::status(ss.str());

   Parameters par = getParameters(0);
}
Parameters ParameterFile::getParameters(int iTime) const {
   if(mParameters.size() == 1) {
      // Assume only one set of parameters for all hours
      std::map<int, Parameters>::const_iterator it = mParameters.begin();
      return it->second;
   }
   else {
      std::map<int, Parameters>::const_iterator it = mParameters.find(iTime);
      if(it == mParameters.end()) {
         std::stringstream ss;
         ss << "Parameter file '" << mFilename << "' does not have values for time " << iTime << ".";
         Util::error(ss.str());
      }
      return it->second;
   }
}

void ParameterFile::setParameters(Parameters iParameters, int iTime) {
   mParameters[iTime] = iParameters;
}

int ParameterFile::getSize() const {
   return mParameters.size();
}
