#include "ParameterFile.h"
#include <fstream>
#include <sstream>
#include "Util.h"

ParameterFile::ParameterFile(std::string iFilename) : 
      mFilename(iFilename) {
   std::ifstream ifs(iFilename.c_str(), std::ifstream::in);
   if(!ifs.good()) {
      Util::error("Parameter file '" + iFilename + "' does not exist");
   }
   while(ifs.good()) {
      char line[10000];
      ifs.getline(line, 10000, '\n');
      if(ifs.good() && line[0] != '#') {
         std::stringstream ss(line);
         // Loop over each value
         std::vector<float> values;
         int time;
         ss >> time;
         while(ss.good()) {
            float value;
            ss >> value;
            values.push_back(value);
         }
         if(values.size() != mNumParameters) {
            std::stringstream ss;
            ss << "Parameter file '" + iFilename + "' is corrupt." << std::endl;
            ss << "The line '" << line << "' does not have " << mNumParameters + 1 << " columns.";
            Util::error(ss.str());
         }
         Parameters parameters(values);
         mParameters[time] = parameters;
      }
   }
   std::stringstream ss;
   ss << "Reading " << mFilename << ". Found " << mParameters.size() << " parameter sets.";
   Util::status(ss.str());
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
