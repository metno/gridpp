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

void writeUsage() {
   std::cout << "Post-processes gridded forecasts" << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  gridpp inputs [options] outputs [options] [-v var [options] [-d downscaler [options] [-p parameters [options]]] [-c calibrator [options] [-p parameters [options]]]*]+" << std::endl;
   std::cout << "        gridpp [--version]" << std::endl;
   std::cout << "        gridpp [--help]" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << "   inputs        One or more input files. Use comma to separate names. When" << std::endl;
   std::cout << "                 string is surrounded by double quotes, wildcards can be used (e.g. *)." << std::endl;
   std::cout << "   outputs       One or more output files. File must already exist. Must" << std::endl;
   std::cout << "                 contain lat/lon information. If multiple input files are" << std::endl;
   std::cout << "                 used, then the same number of outputs must be used." << std::endl;
   std::cout << "   -v var        One of the variables below." << std::endl;
   std::cout << "   -d downscaler One of the downscalers below." << std::endl;
   std::cout << "   -c calibrator One of the calibrators below." << std::endl;
   std::cout << "   -p parameters One of the parameter formats below." << std::endl;
   std::cout << "   options       Options of the form key=value" << std::endl;
   std::cout << "   --version     Print the program's version" << std::endl;
   std::cout << "   --help        Print usage information" << std::endl;
   std::cout << std::endl;
   std::cout << "Notes:" << std::endl;
   std::cout << "   - At least one variable must be specified." << std::endl;
   std::cout << "   - If multiple downscalers are specified for one variable, the last is used." << std::endl;
   std::cout << "   - If the same variable is specified multiple times, the first definition is used." << std::endl;
   std::cout << "   - Multiple identical calibrators are allowed for a single variable." << std::endl;
   std::cout << "   - Only one parameter format can be specified for given a downscaler or calibrator." << std::endl;
   std::cout << std::endl;
   std::cout << "Example:" << std::endl;
   std::cout << "   ./gridpp testing/files/10x10.nc testing/files/10x10_copy.nc -v T -d nearestNeighbour" << std::endl;
   std::cout << std::endl;
   std::cout << "Inputs/Outputs:" << std::endl;
   std::cout << "   I/O types are autodetected, but can be specified using:" << std::endl;
   std::cout << File::getDescriptions();
   std::cout << std::endl;
   std::cout << "Variable options (and default values):" << std::endl;
   std::cout << Util::formatDescription("write=1", "Set to 0 to prevent the variable to be written to output") << std::endl;
   std::cout << Util::formatDescription("units=undef", "Write this to the units attribute in the output.") << std::endl;
   std::cout << Util::formatDescription("standardName=1", "Write this to the standard_name attribute in the output.") << std::endl;
   std::cout << std::endl;
   std::cout << "Downscalers with options (and default values):" << std::endl;
   std::cout << Downscaler::getDescriptions();
   std::cout << std::endl;
   std::cout << "Calibrators with options (and default values):" << std::endl;
   std::cout << Calibrator::getDescriptions();
   std::cout << std::endl;
   std::cout << "Parameter formats with Options (and default values)" << std::endl;
   std::cout << ParameterFile::getDescriptions();
}

int main(int argc, const char *argv[]) {
   double start = Util::clock();

   // Parse command line attributes
   for(int i = 1; i < argc; i++) {
      if(strcmp(argv[i], "--version") == 0) {
         std::cout << "gridpp version " << Util::gridppVersion() << std::endl;
         return 0;
      }
      else if(strcmp(argv[i], "--help") == 0) {
         writeUsage();
         return 0;
      }

   }
   if(argc < 3) {
      writeUsage();
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
   for(int f = 0; f < setup.inputFiles.size(); f++) {
      std::cout << "Input type:  " << setup.inputFiles[f]->name() << std::endl;
      std::cout << "Output type: " << setup.outputFiles[f]->name() << std::endl;
      setup.outputFiles[f]->setTimes(setup.inputFiles[f]->getTimes());
      setup.outputFiles[f]->setReferenceTime(setup.inputFiles[f]->getReferenceTime());

      // Post-process file
      std::vector<Variable> writeVariables;
      for(int v = 0; v < setup.variableConfigurations.size(); v++) {
         double s = Util::clock();
         VariableConfiguration varconf = setup.variableConfigurations[v];
         Variable inputVariable = varconf.inputVariable;
         Variable outputVariable = varconf.outputVariable;

         bool write = 1;
         varconf.outputVariableOptions.getValue("write", write);
         if(write) {
            writeVariables.push_back(outputVariable);
         }
         setup.outputFiles[f]->initNewVariable(outputVariable);

         std::cout << "Processing " << outputVariable.name() << std::endl;

         // Downscale
         std::cout << "   Downscaler " << varconf.downscaler->name() << std::endl;
         varconf.downscaler->downscale(*setup.inputFiles[f], *setup.outputFiles[f]);
         for(int c = 0; c < varconf.calibrators.size(); c++) {
            // Calibrate
            std::cout << "   Calibrator " << varconf.calibrators[c]->name() << std::endl;
            varconf.calibrators[c]->calibrate(*setup.outputFiles[f], varconf.parameterFileCalibrators[c]);
         }
         double e = Util::clock();
         std::cout << "   " << e-s << " seconds" << std::endl;
         std::cout << "Current mem usage input: " << setup.inputFiles[f]->getCacheSize() / 1e6<< std::endl;
         std::cout << "Current mem usage output: " << setup.outputFiles[f]->getCacheSize() / 1e6<< std::endl;
         // setup.inputFile->clear();
      }

      // Write to output
      double s = Util::clock();
      setup.outputFiles[f]->write(writeVariables);
      double e = Util::clock();
      std::cout << "Writing file: " << e-s << " seconds" << std::endl;
      std::cout << "Total time:   " << e-start << " seconds" << std::endl;
      setup.inputFiles[f]->clear();
      setup.outputFiles[f]->clear();
   }
   return 0;
}
