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

   if(iName == "zaga") {
      CalibratorZaga* c = new CalibratorZaga(iVariable, iOptions);

      return c;
   }
   else if(iName == "cloud") {
      CalibratorCloud* c = new CalibratorCloud(iVariable, iOptions);
      return c;
   }
   else if(iName == "accumulate") {
      CalibratorAccumulate* c = new CalibratorAccumulate(iVariable, iOptions);
      return c;
   }
   else if(iName == "gaussian") {
      CalibratorGaussian* c = new CalibratorGaussian(iVariable, iOptions);
      return c;
   }
   else if(iName == "neighbourhood") {
      CalibratorNeighbourhood* c = new CalibratorNeighbourhood(iVariable, iOptions);
      return c;
   }
   else if(iName == "diagnoseDewpoint") {
      CalibratorDiagnoseDewpoint* c = new CalibratorDiagnoseDewpoint(iVariable, iOptions);
      return c;
   }
   else if(iName == "diagnoseRh") {
      CalibratorDiagnoseRh* c = new CalibratorDiagnoseRh(iVariable, iOptions);
      return c;
   }
   else if(iName == "diagnoseWind") {
      CalibratorDiagnoseWind* c = new CalibratorDiagnoseWind(iVariable, iOptions);
      return c;
   }
   else if(iName == "phase") {
      CalibratorPhase* c = new CalibratorPhase(iVariable, iOptions);
      float minPrecip;
      if(iOptions.getValue("minPrecip", minPrecip)) {
         c->setMinPrecip(minPrecip);
      }

      return c;
   }
   else if(iName == "windDirection") {
      CalibratorWindDirection* c = new CalibratorWindDirection(iVariable, iOptions);

      return c;
   }
   else if(iName == "kriging") {
      CalibratorKriging* c = new CalibratorKriging(iVariable, iOptions);

      return c;
   }
   else if(iName == "window") {
      CalibratorWindow* c = new CalibratorWindow(iVariable, iOptions);

      return c;
   }
   else if(iName == "qc") {
      CalibratorQc* c = new CalibratorQc(iVariable, iOptions);

      return c;
   }
   else if(iName == "qnh") {
      CalibratorQnh* c = new CalibratorQnh(iVariable, iOptions);

      return c;
   }
   else if(iName == "qq") {
      CalibratorQq* c = new CalibratorQq(iVariable, iOptions);

      return c;
   }
   else if(iName == "regression") {
      CalibratorRegression* c = new CalibratorRegression(iVariable, iOptions);

      return c;
   }
   else if(iName == "bct") {
      CalibratorBct* c = new CalibratorBct(iVariable, iOptions);

      return c;
   }
   else if(iName == "sort") {
      CalibratorSort* c = new CalibratorSort(iVariable, iOptions);

      return c;
   }
   else if(iName == "altitude") {
      CalibratorAltitude* c = new CalibratorAltitude(iVariable, iOptions);

      return c;
   }
   else if(iName == "mask") {
      CalibratorMask* c = new CalibratorMask(iVariable, iOptions);

      return c;
   }
   else if(iName == "coastal") {
      CalibratorCoastal* c = new CalibratorCoastal(iVariable, iOptions);

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

std::string Calibrator::getDescriptions() {
   std::stringstream ss;
   ss << CalibratorAccumulate::description() << std::endl;
   ss << CalibratorAltitude::description() << std::endl;
   ss << CalibratorBct::description() << std::endl;
   ss << CalibratorCloud::description() << std::endl;
   ss << CalibratorCoastal::description() << std::endl;
   ss << CalibratorDiagnoseDewpoint::description() << std::endl;
   ss << CalibratorDiagnoseRh::description() << std::endl;
   ss << CalibratorDiagnoseWind::description() << std::endl;
   ss << CalibratorGaussian::description() << std::endl;
   ss << CalibratorKriging::description() << std::endl;
   ss << CalibratorMask::description() << std::endl;
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
