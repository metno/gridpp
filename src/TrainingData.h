#ifndef TRAININGDATA_H
#define TRAININGDATA_H
#include <string>
#include <vector>
#include "Location.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <map>

typedef std::vector<float> Ens;
typedef std::pair<float,Ens> ObsEns;
class TrainingData {
   public:
      TrainingData(std::string iFilename);
      // std::vector<float>               getObs(const Location& iLocation);
      // std::vector<std::vector<float> > getEns(const Location& iLocation);
      std::vector<ObsEns> getData(int iOffset) const;
      std::vector<int> getOffsets() const;
      std::string getFilename() const;
   private:
      std::string mFilename;
      static float getRand();
      std::map<int, std::vector<ObsEns> > mData; // offset, data
      static boost::variate_generator<boost::mt19937, boost::uniform_01<> > mRand;
};
#endif
