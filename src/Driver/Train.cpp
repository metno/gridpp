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
   std::cout << "Trains a calibration method in gridpp" << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  gridpp_train dataFile -p parameters [options] -c calibrator [options] -v variable" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << Util::formatDescription("dataFile=required","Input text file with obs/forecast data. Must be in the following format:") << std::endl;
   std::cout << Util::formatDescription("","<time> <lat> <lon> <elev> <obs> <ens1> <ens2> ... <ensN>") << std::endl;
   std::cout << Util::formatDescription("-p parameters","Put parameters into this output format.") << std::endl;
   std::cout << Util::formatDescription("-c calibrator","Train this calibration method") << std::endl;
   std::cout << Util::formatDescription("-v variable","Train this variable") << std::endl;
   std::cout << std::endl;
   std::cout << "Example:" << std::endl;
   std::cout << "   ./gridpp_train testing/files/training.txt -p text file=parameters.txt -c regression -v T" << std::endl;
   std::cout << std::endl;
   std::cout << "Variables:" << std::endl;
   std::cout << Variable::getDescriptions();
   std::cout << std::endl;
   std::cout << "Calibrators with options (and default values):" << std::endl;
   std::cout << Calibrator::getDescriptions();
   std::cout << std::endl;
   std::cout << "Parameter formats with Options (and default values)" << std::endl;
   std::cout << ParameterFile::getDescriptions();
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
      setup.output->setParameters(par, offset, Location(Util::MV, Util::MV, Util::MV));
   }

   setup.output->write();
}
