#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
#include "../Options.h"
#include "../Setup.h"
#include "../TrainingData.h"
void writeUsage() {
   std::cout << "Trains a method in gridpp" << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  train data method options" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << Util::formatDescription("data=required","Input text file with obs/forecast data. Must be in the following format:") << std::endl;
   std::cout << Util::formatDescription("","<time> <lat> <lon> <elev> <obs> <fcst>") << std::endl;
   // std::cout << Util::formatDescription("output=undef","Bias file, suitable for kriging.") << std::endl;
   std::cout << Util::formatDescription("method=undef","Current Kalman coefficients. If unspecified, initialize. Must have the format:") << std::endl;
   std::cout << Util::formatDescription("","<time> <lat> <lon> <elev> <pVarV> <kalmanGainVar> <varV> <p> <kalmanGain> <error> <lastError> <biasEstimate>") << std::endl;
   std::cout << std::endl;
}
int main(int argc, const char *argv[]) {
   // if(argc <= 1) {
   //    writeUsage();
   //    return 0;
   // }

   std::string dataFilename = "";
   std::string outputFilename = "test.txt";
   std::string newFilename = "";
   if(argc > 1)
      dataFilename = argv[1];
   /*
   Options options;
   for(int i = 1; i < argc; i++) {
      options.addOptions(std::string(argv[i]));
   }
   if(!options.getValue("data", dataFilename)) {
      Util::error("kalmanFilter: 'data' required.");
   }
   */
   Options dataOptions("file=testing/files/obsfcst.txt spatial=1");
   Options outputOptions("file=newParameters.txt spatial=0");
   Options schemeOptions("fractThreshold=1");
   Variable::Type variable = Variable::T;

   ParameterFile* dataFile   = ParameterFile::getScheme("text", dataOptions);
   ParameterFile* outputFile = ParameterFile::getScheme("text", outputOptions);

   // Scheme* scheme = new CalibratorZaga(outputFile, variable, schemeOptions);
   CalibratorZaga* scheme = new CalibratorZaga(outputFile, variable, schemeOptions);
   TrainingData trainingData(dataFilename);
   for(int i = 6; i < 12; i += 6) {
      Parameters par = scheme->train(trainingData, i);
      std::cout << "Parameters time " << i << ":";
      for(int k = 0; k < par.size(); k++) {
         std::cout << " " << par[k];
      }
      std::cout << std::endl;
      outputFile->setParameters(par, i/6-1, Location(Util::MV, Util::MV, Util::MV));
   }

   outputFile->write();
}
