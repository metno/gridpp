#include "TrainingData.h"
#include <fstream>
#include <sstream>
boost::variate_generator<boost::mt19937, boost::uniform_01<> > TrainingData::mRand(boost::mt19937(0), boost::uniform_01<>());

TrainingData::TrainingData(std::string iFilename) : mFilename(iFilename) {
   if(mFilename != "") {
      std::ifstream ifs(mFilename.c_str(), std::ifstream::in);
      if(!ifs.good()) {
         Util::error("File '" + mFilename + "' does not exist");
      }
      int counter = 0;
      char line[10000];
      while(ifs.good() && mData.size() < 10000) {
         ifs.getline(line, 10000, '\n');
         if(ifs.good() && line[0] != '#') {
            std::stringstream ss(line);
            // Loop over each value
            std::vector<float> ens;
            int id, date, init, offset;
            float lat, lon, elev, obs;
            // ss >> id >> lat >> lon >> elev >> date >> init >> offset >> obs;
            ss >> offset >> lat >> lon >> elev >> obs;
            bool isValid = Util::isValid(obs);
            while(ss.good()) {
               float value;
               bool status  = ss >> value;
               if(!Util::isValid(value)) {
                  isValid = false;
               }
               if(status) {
                  ens.push_back(value);
               }
            }
            if(isValid) {
               ObsEns obsEns(obs, ens);
               mData[offset].push_back(obsEns);
               counter++;
            }
         }
      }
   }
   else {
      for(int t = 0; t < 1000; t++) {
         float obs = 3.2;
         std::vector<float> ens(5,0);
         for(int j = 0; j < 5; j++) {
            ens[j] = getRand();
         }
         ObsEns obsEns(obs, ens);
         mData[0].push_back(obsEns);
      }
   }
}

std::vector<ObsEns> TrainingData::getData(int iOffset) const {
   std::map<int, std::vector<ObsEns> >::const_iterator it = mData.find(iOffset);
   if(it == mData.end()) {
      std::vector<ObsEns> empty;
      return empty;
   }
   else {
      return it->second;
   }
}

float TrainingData::getRand() {
   float value = mRand();
   if(value < 0.1)
      return 0;
   else
      return (value-0.1)*5;
}

std::vector<int> TrainingData::getOffsets() const {
   std::vector<int> offsets;
   std::map<int, std::vector<ObsEns> >::const_iterator it;
   for(it = mData.begin(); it != mData.end(); it++) {
      offsets.push_back(it->first);
   }
   return offsets;
}

std::string TrainingData::getFilename() const {
   return mFilename;
}
