#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
#include "../Options.h"
#include "../SetupTrain.h"
#include "../TrainingData.h"
void writeUsage() {
   std::cout << "Trains a method in gridpp" << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  gridpp_train dataFile -o output [options] -c method [options] -v variable" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << Util::formatDescription("dataFile=required","Input text file with obs/forecast data. Must be in the following format:") << std::endl;
   std::cout << Util::formatDescription("","<time> <lat> <lon> <elev> <obs> <fcst>") << std::endl;
   std::cout << Util::formatDescription("-o output","Put parameters into this output file") << std::endl;
   std::cout << Util::formatDescription("-m method","Train this calibration method") << std::endl;
   std::cout << Util::formatDescription("-v variable","Train this variable") << std::endl;
   std::cout << std::endl;
}
int main(int argc, const char *argv[]) {
   if(argc <= 1) {
      writeUsage();
      return 0;
   }
   // Retrieve setup
   std::vector<std::string> args;
   for(int i = 1; i < argc; i++) {
      args.push_back(std::string(argv[i]));
   }
   SetupTrain setup(args);

   Variable::Type variable = Variable::T;
   std::vector<int> offsets = setup.trainingData->getOffsets();

   for(int i = 0; i < offsets.size(); i++) {
      int offset = offsets[i];
      double timeStart = Util::clock();
      Parameters par = setup.method->train(*setup.trainingData, offset);
      double timeEnd = Util::clock();
      std::cout << "Time: " << timeEnd - timeStart << std::endl;
      std::cout << "Calculating parameters for offset=" << i << ": ";
      for(int k = 0; k < par.size(); k++) {
         std::cout << " " << par[k];
      }
      std::cout << std::endl;
      setup.output->setParameters(par, i/6-1, Location(Util::MV, Util::MV, Util::MV));
   }

   setup.output->write();
}
