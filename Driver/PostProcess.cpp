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
      std::cout << "usage:  postprocess.exe input output [-v var [-d downscaler [options]*]] [-c calibrator [options]*]]*]+" << std::endl;
      std::cout << std::endl;
      std::cout << "Arguments:" << std::endl;
      std::cout << "   input         Netcdf file with AROME data." << std::endl;
      std::cout << "   output        Netcdf file with AROME data." << std::endl;
      std::cout << "   -v var        Variable." << std::endl;
      std::cout << "   -d downscaler One of the downscalers below." << std::endl;
      std::cout << "   -c calibrator One of the calibrators below." << std::endl;
      std::cout << "   options       Options of the form key=value" << std::endl;

      std::cout << std::endl;
      std::cout << "Variables:" << std::endl;
      std::cout << Variable::description();
      std::cout << std::endl;
      std::cout << "Downscalers with options (and default values):" << std::endl;
      std::cout << DownscalerNearestNeighbour::description();
      std::cout << DownscalerGradient::description();
      std::cout << DownscalerSmart::description();
      std::cout << std::endl;
      std::cout << "Calibrators with options (and default values):" << std::endl;
      std::cout << CalibratorZaga::description();
      std::cout << CalibratorCloud::description();
      std::cout << CalibratorAccumulate::description();
      return 1;
   }
   Util::setShowError(true);
   Util::setShowWarning(true);
   Util::setShowStatus(true);

   // Retrieve setup
   Setup setup;
   std::vector<std::string> args;
   for(int i = 1; i < argc; i++) {
      args.push_back(argv[i]);
   }
   bool success = Setup::getSetup(args, setup);
   if(!success) {
      Util::error("Could not configure setup");
   }

   // Post-process file
   std::vector<Variable::Type> variables;
   for(int v = 0; v < setup.variableConfigurations.size(); v++) {
      double s = Util::clock();
      VariableConfiguration varconf = setup.variableConfigurations[v];
      Variable::Type variable = varconf.variable;
      variables.push_back(variable);
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
   }

   // Write to output
   double s = Util::clock();
   setup.outputFile->write(variables);
   double e = Util::clock();
   std::cout << "Writing file: " << e-s << " seconds" << std::endl;
   std::cout << "Total time:   " << e-start << " seconds" << std::endl;
}
