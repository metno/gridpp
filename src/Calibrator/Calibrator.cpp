#include "Calibrator.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../Options.h"
#include "../ParameterFile/ParameterFile.h"
#include "../File/File.h"

Calibrator::Calibrator(const Variable& iVariable, const Options& iOptions) : Scheme(iOptions),
      mVariable(iVariable),
      mOptions(iOptions) {

}
Calibrator* Calibrator::getScheme(std::string iName, Variable iVariable, const Options& iOptions) {
   Calibrator* c;

   if(iName == "zaga") {
      c = new CalibratorZaga(iVariable, iOptions);
   }
   else if(iName == "cloud") {
      c = new CalibratorCloud(iVariable, iOptions);
   }
   else if(iName == "accumulate") {
      c = new CalibratorAccumulate(iVariable, iOptions);
   }
   else if(iName == "deaccumulate") {
      c = new CalibratorDeaccumulate(iVariable, iOptions);
   }
   else if(iName == "gaussian") {
      c = new CalibratorGaussian(iVariable, iOptions);
   }
   else if(iName == "neighbourhood") {
      c = new CalibratorNeighbourhood(iVariable, iOptions);
   }
   else if(iName == "diagnoseHumidity") {
      c = new CalibratorDiagnoseHumidity(iVariable, iOptions);
   }
   else if(iName == "diagnoseWind") {
      c = new CalibratorDiagnoseWind(iVariable, iOptions);
   }
   else if(iName == "oi") {
      c = new CalibratorOi(iVariable, iOptions);
   }
   else if(iName == "override") {
      c = new CalibratorOverride(iVariable, iOptions);
   }
   else if(iName == "phase") {
      c = new CalibratorPhase(iVariable, iOptions);
   }
   else if(iName == "windDirection") {
      c = new CalibratorWindDirection(iVariable, iOptions);
   }
   else if(iName == "kriging") {
      c = new CalibratorKriging(iVariable, iOptions);
   }
   else if(iName == "window") {
      c = new CalibratorWindow(iVariable, iOptions);
   }
   else if(iName == "qc") {
      c = new CalibratorQc(iVariable, iOptions);
   }
   else if(iName == "qnh") {
      c = new CalibratorQnh(iVariable, iOptions);
   }
   else if(iName == "qq") {
      c = new CalibratorQq(iVariable, iOptions);
   }
   else if(iName == "regression") {
      c = new CalibratorRegression(iVariable, iOptions);
   }
   else if(iName == "bct") {
      c = new CalibratorBct(iVariable, iOptions);
   }
   else if(iName == "sort") {
      c = new CalibratorSort(iVariable, iOptions);
   }
   else if(iName == "altitude") {
      c = new CalibratorAltitude(iVariable, iOptions);
   }
   else if(iName == "mask") {
      c = new CalibratorMask(iVariable, iOptions);
   }
   else if(iName == "coastal") {
      c = new CalibratorCoastal(iVariable, iOptions);
   }
   else {
      Util::error("Could not instantiate calibrator with name '" + iName + "'");
      return NULL;
   }
   return c;
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
   // Arrange data
   std::vector<ObsEns> data;
   for(int d = 0; d < iData.size(); d++){
      // Downscaling (currently nearest neighbour)
      float obs = (*(iData[d].first))(iIobs,iJobs,0);
      Ens ens   = (*(iData[d].second))(iIens, iJens);
      ObsEns obsens(obs, ens);
      data.push_back(obsens);
   }

   Parameters par = train(data);
   return par;
}

Options Calibrator::getOptions() const {
   return mOptions;
}

std::string Calibrator::getDescriptions(bool full) {
   std::stringstream ss;
   ss << CalibratorAccumulate::description(full) << std::endl;
   ss << CalibratorAltitude::description(full) << std::endl;
   ss << CalibratorBct::description(full) << std::endl;
   ss << CalibratorCloud::description(full) << std::endl;
   ss << CalibratorCoastal::description(full) << std::endl;
   ss << CalibratorDeaccumulate::description(full) << std::endl;
   ss << CalibratorDiagnoseHumidity::description(full) << std::endl;
   ss << CalibratorDiagnoseWind::description(full) << std::endl;
   ss << CalibratorGaussian::description(full) << std::endl;
   ss << CalibratorKriging::description(full) << std::endl;
   ss << CalibratorMask::description(full) << std::endl;
   ss << CalibratorNeighbourhood::description(full) << std::endl;
   ss << CalibratorOi::description(full) << std::endl;
   ss << CalibratorOverride::description(full) << std::endl;
   ss << CalibratorPhase::description(full) << std::endl;
   ss << CalibratorQc::description(full) << std::endl;
   ss << CalibratorQnh::description(full) << std::endl;
   ss << CalibratorQq::description(full) << std::endl;
   ss << CalibratorRegression::description(full) << std::endl;
   ss << CalibratorSort::description(full) << std::endl;
   ss << CalibratorWindDirection::description(full) << std::endl;
   ss << CalibratorZaga::description(full) << std::endl;
   return ss.str();
}
