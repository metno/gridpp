#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
#include "../Options.h"
#include "../Setup.h"

int main(int argc, const char *argv[]) {
   double start = Util::clock();

   // Parse command line attributes
   if(argc < 3) {
      std::cout << "Post-processes gridded forecasts" << std::endl;
      std::cout << std::endl;
      std::cout << "usage:  gridpp inFilename [options]* outFilename [options]* [-v var [options]* [-d downscaler [options]*]] [-c calibrator [options]*]]*]+" << std::endl;
      std::cout << std::endl;
      std::cout << "Arguments:" << std::endl;
      std::cout << "   input         Input file with NetCDF data." << std::endl;
      std::cout << "   output        Output file with NetCDF data. File must already exist. Must" << std::endl;
      std::cout << "                 contain lat/lon information." << std::endl;
      std::cout << "   -v var        One of the variables below." << std::endl;
      std::cout << "   -d downscaler One of the downscalers below." << std::endl;
      std::cout << "   -c calibrator One of the calibrators below." << std::endl;
      std::cout << "   options       Options of the form key=value" << std::endl;
      std::cout << std::endl;
      std::cout << "Notes:" << std::endl;
      std::cout << "   - At least one variable must be specified." << std::endl;
      std::cout << "   - If multiple downscalers are specified for one variable, the last is used." << std::endl;
      std::cout << "   - If the same variable is specified multiple times, the first definition is used." << std::endl;
      std::cout << "   - Multiple identical calibrators are allowed for a single variable." << std::endl;

      std::cout << std::endl;
      std::cout << "Inputs/Outputs:" << std::endl;
      std::cout << "   I/O types are autodetected, but can be specified using:" << std::endl;
      std::cout << FileArome::description();
      std::cout << FileEc::description();
      std::cout << FilePoint::description();
      std::cout << std::endl;
      std::cout << "Variables:" << std::endl;
      std::cout << Variable::description();
      std::cout << std::endl;
      std::cout << "Variable options (and default values):" << std::endl;
      std::cout << "   write=1                      Set to 0 to prevent the variable to be written to output" << std::endl;
      std::cout << std::endl;
      std::cout << "Downscalers with options (and default values):" << std::endl;
      std::cout << DownscalerNearestNeighbour::description();
      std::cout << DownscalerGradient::description();
      std::cout << DownscalerSmart::description();
      std::cout << DownscalerPressure::description();
      std::cout << DownscalerBypass::description();
      std::cout << std::endl;
      std::cout << "Calibrators with options (and default values):" << std::endl;
      std::cout << CalibratorZaga::description();
      std::cout << CalibratorCloud::description();
      std::cout << CalibratorQnh::description();
      std::cout << CalibratorAccumulate::description();
      std::cout << CalibratorWindDirection::description();
      std::cout << CalibratorSmooth::description();
      std::cout << CalibratorPhase::description();
      return 0;
   }
   Util::setShowError(true);
   Util::setShowWarning(true);
   Util::setShowStatus(false);

   // Retrieve setup
   std::vector<std::string> args;
   for(int i = 1; i < argc; i++) {
      args.push_back(std::string(argv[i]));
   }
   Setup setup(args);
   std::cout << "Input type:  " << setup.inputFile->name() << std::endl;
   std::cout << "Output type: " << setup.outputFile->name() << std::endl;
   setup.outputFile->setTimes(setup.inputFile->getTimes());
   setup.outputFile->setReferenceTime(setup.inputFile->getReferenceTime());

   // Post-process file
   std::vector<Variable::Type> writeVariables;
   for(int v = 0; v < setup.variableConfigurations.size(); v++) {
      double s = Util::clock();
      VariableConfiguration varconf = setup.variableConfigurations[v];
      Variable::Type variable = varconf.variable;

      bool write = 1;
      varconf.variableOptions.getValue("write", write);
      if(write) {
         writeVariables.push_back(variable);
      }
      setup.outputFile->initNewVariable(variable);

      std::cout << "Processing " << Variable::getTypeName(variable) << std::endl;

      // Downscale
      std::cout << "   Downscaler " << varconf.downscaler->name() << std::endl;
      varconf.downscaler->downscale(*setup.inputFile, *setup.outputFile);
      for(int c = 0; c < varconf.calibrators.size(); c++) {
         // Calibrate
         std::cout << "   Calibrator " << varconf.calibrators[c]->name() << std::endl;
         varconf.calibrators[c]->calibrate(*setup.outputFile);
      }
      double e = Util::clock();
      std::cout << "   " << e-s << " seconds" << std::endl;
      std::cout << "Current mem usage input: " << setup.inputFile->getCacheSize() / 1e6<< std::endl;
      std::cout << "Current mem usage output: " << setup.outputFile->getCacheSize() / 1e6<< std::endl;
      // setup.inputFile->clear();
   }

   // Write to output
   double s = Util::clock();
   setup.outputFile->write(writeVariables);
   double e = Util::clock();
   std::cout << "Writing file: " << e-s << " seconds" << std::endl;
   std::cout << "Total time:   " << e-start << " seconds" << std::endl;

   return 0;
}
