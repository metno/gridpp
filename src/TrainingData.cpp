#include "TrainingData.h"
#include <fstream>
#include <sstream>
boost::variate_generator<boost::mt19937, boost::uniform_01<> > TrainingData::mRand(boost::mt19937(0), boost::uniform_01<>());

TrainingData::TrainingData(std::string iFilename) : mFilename(iFilename) {

}

std::vector<ObsEns> TrainingData::getData(int iOffset) const {
   std::vector<ObsEns> data;
   if(mFilename != "") {
      std::cout << "Reading " << mFilename << std::endl;
      std::ifstream ifs(mFilename.c_str(), std::ifstream::in);
      if(!ifs.good()) {
         Util::error("File '" + mFilename + "' does not exist");
      }
      int counter = 0;
      char line[10000];
      ifs.getline(line, 10000, '\n');
      while(ifs.good() && data.size() < 10000) {
         ifs.getline(line, 10000, '\n');
         if(ifs.good() && line[0] != '#') {
            std::stringstream ss(line);
            // Loop over each value
            std::vector<float> ens;
            int id, date, init, offset;
            float lat, lon, elev, obs;
            ss >> id >> lat >> lon >> elev >> date >> init >> offset >> obs;
            if(offset == iOffset) {
               while(ss.good()) {
                  float value;
                  bool status  = ss >> value;
                  if(status) {
                     ens.push_back(value);
                  }
               }
               ObsEns obsEns(obs, ens);
               data.push_back(obsEns);
               counter++;
            }
         }
      }
      std::cout << "Found " << data.size() << " rows" << std::endl;
   }
   else {
      for(int t = 0; t < 1000; t++) {
         float obs = 3.2;
         std::vector<float> ens(5,0);
         for(int j = 0; j < 5; j++) {
            ens[j] = getRand();
         }
         ObsEns obsEns(obs, ens);
         data.push_back(obsEns);
      }
   }
   return data;
}

float TrainingData::getRand() {
   float value = mRand();
   if(value < 0.1)
      return 0;
   else
      return (value-0.1)*5;
}
