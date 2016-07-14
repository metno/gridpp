#include "Calibrator.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../Options.h"
#include "../ParameterFile/ParameterFile.h"
#include "../File/File.h"

Calibrator::Calibrator(const Options& iOptions) : Scheme(iOptions) {

}
Calibrator* Calibrator::getScheme(std::string iName, const Options& iOptions) {

   if(iName == "zaga") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'zaga' needs variable");
      }
      CalibratorZaga* c = new CalibratorZaga(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "cloud") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'cloud' needs variable");
      }
      CalibratorCloud* c = new CalibratorCloud(Variable::getType(variable), iOptions);
      return c;
   }
   else if(iName == "accumulate") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'accumulate' needs variable");
      }
      CalibratorAccumulate* c = new CalibratorAccumulate(Variable::getType(variable), iOptions);
      return c;
   }
   else if(iName == "gaussian") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'gaussian' needs variable");
      }
      CalibratorGaussian* c = new CalibratorGaussian(Variable::getType(variable), iOptions);
      return c;
   }
   else if(iName == "neighbourhood") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'neighbourhood' needs variable");
      }
      CalibratorNeighbourhood* c = new CalibratorNeighbourhood(Variable::getType(variable), iOptions);
      return c;
   }
   else if(iName == "phase") {
      CalibratorPhase* c = new CalibratorPhase(iOptions);
      float minPrecip;
      if(iOptions.getValue("minPrecip", minPrecip)) {
         c->setMinPrecip(minPrecip);
      }
      bool useWetbulb;
      if(iOptions.getValue("useWetbulb", useWetbulb)) {
         c->setUseWetbulb(useWetbulb);
      }

      return c;
   }
   else if(iName == "windDirection") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'windDirection' needs variable");
      }
      CalibratorWindDirection* c = new CalibratorWindDirection(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "kriging") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'kriging' needs variable");
      }
      CalibratorKriging* c = new CalibratorKriging(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "diagnose") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'diagnose' needs variable");
      }
      CalibratorDiagnose* c = new CalibratorDiagnose(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "window") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'window' needs variable");
      }
      CalibratorWindow* c = new CalibratorWindow(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "qc") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'qc' needs variable");
      }
      CalibratorQc* c = new CalibratorQc(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "qnh") {
      CalibratorQnh* c = new CalibratorQnh(iOptions);

      return c;
   }
   else if(iName == "qq") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'regression' needs variable");
      }
      CalibratorQq* c = new CalibratorQq(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "regression") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'regression' needs variable");
      }
      CalibratorRegression* c = new CalibratorRegression(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "bct") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'bct' needs variable");
      }
      CalibratorBct* c = new CalibratorBct(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "sort") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'sort' needs variable");
      }
      CalibratorSort* c = new CalibratorSort(Variable::getType(variable), iOptions);

      return c;
   }
   else if(iName == "altitude") {
      CalibratorAltitude* c = new CalibratorAltitude(iOptions);

      return c;
   }
   else if(iName == "coastal") {
      std::string variable;
      if(!iOptions.getValue("variable", variable)) {
         Util::error("Calibrator 'sort' needs variable");
      }
      CalibratorCoastal* c = new CalibratorCoastal(Variable::getType(variable), iOptions);

      return c;
   }
   else {
      Util::error("Could not instantiate calibrator with name '" + iName + "'");
      return NULL;
   }
}
bool Calibrator::calibrate(File& iFile, const ParameterFile* iParameterFile) const {
   if(requiresParameterFile() && iParameterFile == NULL) {
      std::stringstream ss;
      ss << "Calibrator '" << name() << "' requires a parameter file";
      Util::error(ss.str());
   }
   return calibrateCore(iFile, iParameterFile);
}

void Calibrator::shuffle(const std::vector<float>& iBefore, std::vector<float>& iAfter) {
   if(iBefore.size() != iAfter.size()) {
      return;
   }
   if(iBefore.size() == 0)
      return;

   int N = iBefore.size();
   std::vector<std::pair<float,int> > pairs(N);
   for(int e = 0; e < N; e++) {
      if(!Util::isValid(iBefore[e]) || !Util::isValid(iAfter[e])) {
         return;
      }
      pairs[e].first = iBefore[e];
      pairs[e].second = e;
   }
   std::vector<float> afterCopy = iAfter;
   // Sort values so that the rank of a member is the same before and after calibration
   std::sort(pairs.begin(), pairs.end(), Util::sort_pair_first<float,int>());
   std::sort(afterCopy.begin(), afterCopy.end());
   for(int e = 0; e < N; e++) {
      int ei = pairs[e].second;
      float valueCal = afterCopy[e];
      iAfter[ei] = valueCal;
   }
}

Parameters Calibrator::train(const std::vector<ObsEns>& iData) const {
   Util::error("Cannot train method. Not yet implemented.");
   return Parameters();
}

Parameters Calibrator::train(const std::vector<ObsEnsField>& iData, const Grid& iObsGrid, const Grid& iEnsGrid, int iIobs, int iJobs, int iIens, int iJens) const {
   Util::error("Cannot train method. Not yet implemented.");
   // Arrange data
   std::vector<ObsEns> data;
   for(int d = 0; d < iData.size(); d++){
      // Downscaling (currently nearest neighbour)
      int iIobs = 0;
      int iJobs = 0;
      abort();
      float obs = (*(iData[d].first))(iIobs,iJobs,0);
      Ens ens   = (*(iData[d].second))(iIens, iJens);
      ObsEns obsens(obs, ens);
      data.push_back(obsens);
   }

   Parameters par = train(data);
   return par;
}

std::string Calibrator::getDescriptions() {
   std::stringstream ss;
   ss << CalibratorAccumulate::description() << std::endl;
   ss << CalibratorAltitude::description() << std::endl;
   ss << CalibratorBct::description() << std::endl;
   ss << CalibratorCloud::description() << std::endl;
   ss << CalibratorCoastal::description() << std::endl;
   ss << CalibratorDiagnose::description() << std::endl;
   ss << CalibratorGaussian::description() << std::endl;
   ss << CalibratorKriging::description() << std::endl;
   ss << CalibratorNeighbourhood::description() << std::endl;
   ss << CalibratorPhase::description() << std::endl;
   ss << CalibratorQc::description() << std::endl;
   ss << CalibratorQnh::description() << std::endl;
   ss << CalibratorQq::description() << std::endl;
   ss << CalibratorRegression::description() << std::endl;
   ss << CalibratorSort::description() << std::endl;
   ss << CalibratorWindDirection::description() << std::endl;
   ss << CalibratorZaga::description() << std::endl;
   return ss.str();
}
