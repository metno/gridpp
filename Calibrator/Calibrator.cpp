#include "Calibrator.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"

Calibrator::Calibrator(const ParameterFile& iParameterFile):
      mParameterFile(iParameterFile) {
}
void Calibrator::calibrate(DataFile& iFile) const {
   calibrateCore(iFile);
}

void Calibrator::shuffle(const std::vector<float>& iBefore, std::vector<float>& iAfter) {
   assert(iBefore.size() == iAfter.size());
   int N = iBefore.size();
   std::vector<std::pair<float,int> > pairs(N);
   bool isValid = true;
   for(int e = 0; e < N; e++) {
      if(!Util::isValid(iBefore[e])) {
         isValid = false;
         break;
      }
      pairs[e].first = iBefore[e];
      pairs[e].second = e;
   }
   if(isValid) {
      std::vector<float> afterCopy = iAfter;
      // Sort values so that the rank of a member is the same before and after calibration
      std::sort(pairs.begin(), pairs.end(), sort_pair_first<float,int>());
      for(int e = 0; e < N; e++) {
         int ei = pairs[e].second;
         float valueCal = afterCopy[e];
         iAfter[ei] = valueCal;
      }
   }
}

float Calibrator::logit(float p) {
   return log(p/(1-p));
}
float Calibrator::invLogit(float x) {
   return exp(x)/(exp(x)+1);
}
