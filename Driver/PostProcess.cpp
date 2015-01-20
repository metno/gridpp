#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
#include "../Options.h"
struct VariableConfiguration {
   Variable::Type variable;
   Downscaler* downscaler;
   std::vector<Calibrator*> calibrators;
};
struct Setup {
   File* inputFile;
   File* outputFile;
   std::vector<VariableConfiguration> variableConfigurations;
};
bool getSetup(int argc, const char *argv[], Setup& iSetup);

int main(int argc, const char *argv[]) {
   double start = Util::clock();

   // Parse command line attributes
   if(argc < 3) {
      std::cout << "Post-processes gridded forecasts" << std::endl;
      std::cout << std::endl;
      std::cout << "usage:  postprocess.exe input output [-v var [-d downscaler [-dp dParFile]] [-c calibrator [-cp cParFile]]*]+" << std::endl;
      std::cout << std::endl;
      std::cout << "Arguments:" << std::endl;
      std::cout << "   input         Netcdf file with AROME data." << std::endl;
      std::cout << "   output        Netcdf file with AROME data." << std::endl;
      std::cout << "   parameters    Text file with parameters." << std::endl;
      std::cout << "   -v            Verbose. Show status messages." << std::endl;
      return 1;
   }
   Util::setShowError(true);
   Util::setShowWarning(true);
   Util::setShowStatus(false);

   // Retrieve setup
   Setup setup;
   bool success = getSetup(argc, argv, setup);//ifile, ofile, variables, downscalers, calibrators);
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

      // Downscale
      varconf.downscaler->downscale(*setup.inputFile, *setup.outputFile);
      for(int c = 0; c < varconf.calibrators.size(); c++) {
         // Calibrate
         varconf.calibrators[c]->calibrate(*setup.outputFile);
      }
      double e = Util::clock();
      std::cout << "Processing " << Variable::getTypeName(variable) << ": " << e-s << " seconds" << std::endl;
   }

   // Write to output
   double s = Util::clock();
   setup.outputFile->write(variables);
   double e = Util::clock();
   std::cout << "Writing file: " << e-s << " seconds" << std::endl;
   std::cout << "Total time:   " << e-start << " seconds" << std::endl;
}

bool getSetup(std::vector<std::string> argv, Setup& iSetup) {
   // Implement a finite state machine
   enum State {START = 0, VAR = 1, NEWVAR = 2, DOWN = 10, DOWNOPT = 15, CAL = 20, NEWCAL = 22, CALOPT = 25, END = 30, ERROR = 40};
   State state = START;
   Variable::Type variable;
   Options dOptions;
   Options cOptions;
   int index = 0;
   std::string downscaler = "nearestneighbour";
   std::string calibrator = "";
   std::vector<Calibrator*> calibrators;
   while(true) {
      if(state == START) {
         if(argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else {
            state = ERROR;
         }
      }
      else if(state == VAR) {
         variable = Variable::getType(argv[index]);
         index++;
         if(argv.size() == index) {
            state = END;
         }
         else if(argv[index] == "-d") {
            state = DOWN;
         }
         else if(argv[index] == "-c") {
            state = CAL;
         }
         else {
            state = ERROR;
         }
         index++;
      }
      else if(state = NEWVAR) {
         Downscaler* d = Downscaler::getScheme(downscaler, variable, dOptions);
         VariableConfiguration varconf;
         varconf.variable = variable;
         varconf.downscaler = d;
         varconf.calibrators = calibrators;

         // Reset to defaults
         downscaler = "nearestneighbour";
         dOptions.clear();
         calibrators.clear();

         if(argv.size() < index) {
            state = END;
         }
         else {
            state = VAR;
         }
      }
      else if(state == DOWN) {
         downscaler = argv[index];
         index++;
         if(argv[index] == "-c") {
            state = CAL;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
            index++;
         }
         else if(argv.size() < index) {
            state = NEWVAR;
         }
         else {
            state = DOWNOPT;
         }
      }
      else if(state == DOWNOPT) {
         if(argv[index] == "-c") {
            state = CAL;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
            index++;
         }
         else if(argv.size() < index) {
            state = NEWVAR;
         }
         else {
            // Process downscaler options
            if(argv.size() < index+1) {
               dOptions.addOption(argv[index], argv[index+1]);
               index++;
               index++;
            }
            else {
               state = ERROR;
            }
         }
      }
      else if(state == CAL) {
         calibrator = argv[index];
         index++;
         if(argv[index] == "-v") {
            state = NEWVAR;
            index++;
         }
         else if(argv[index] == "-c") {
            state = NEWCAL;
            index++;
         }
         else {
            state = CALOPT;
         }
      }
      else if(state == CALOPT) {
         if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] ==  "-c") {
            state = NEWCAL;
         }
         else {
            // Process calibrator options
            if(argv.size() < index+1) {
               cOptions.addOption(argv[index], argv[index+1]);
               index++;
               index++;
            }
            else {
               state = ERROR;
            }
         }
      }
      else if(state == NEWCAL) {
         Calibrator* c = Calibrator::getScheme(calibrator, variable, cOptions);
         calibrators.push_back(c);

         // Reset
         calibrator = "";
         cOptions.clear();
      }
      else if(state == END) {
         break;
      }
      else if(state == ERROR) {
         Util::error("error in FSM");
      }
   }
}
