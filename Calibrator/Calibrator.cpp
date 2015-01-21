#include "Calibrator.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../Options.h"

Calibrator::Calibrator() {

}
Calibrator* Calibrator::getScheme(std::string iType, Options& iOptions) {

   if(iType == "zaga") {
      std::string parFilename;
      if(!iOptions.getValue("parameters", parFilename)) {
         Util::error("Calibrator 'zaga' needs parameters");
      }

      ParameterFile* parFile = new ParameterFile(parFilename);
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'zaga' needs variable");
      }
      CalibratorZaga* c = new CalibratorZaga(parFile, Variable::getType(variable));

      // Optional settings
      float fracThreshold;
      if(iOptions.getValue("fracThreshold", fracThreshold)) {
         c->setFracThreshold(fracThreshold);
      }
      return c;
   }
   else if(iType == "cloud") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'cloud' needs variable");
      }
      CalibratorCloud* c = new CalibratorCloud(Variable::Precip, Variable::getType(variable));
      return c;
   }
   else if(iType == "accumulate") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'accumulate' needs variable");
      }
      CalibratorAccumulate* c = new CalibratorAccumulate(Variable::getType(variable));
      return c;
   }
   else {
      Util::error("Could not instantiate calibrator of type '" + iType + "'");
      return NULL;
   }
}
void Calibrator::calibrate(File& iFile) const {
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
