#include "Setup.h"
#include "File/File.h"
#include "Calibrator/Calibrator.h"
#include "Downscaler/Downscaler.h"

bool Setup::getSetup(std::vector<std::string> argv, Setup& iSetup) {
   iSetup.inputFile  = File::getScheme(argv[0]);
   iSetup.outputFile = File::getScheme(argv[1]);

   // Implement a finite state machine
   enum State {START = 0, VAR = 1, NEWVAR = 2, DOWN = 10, DOWNOPT = 15, CAL = 20, NEWCAL = 22, CALOPT = 25, END = 30, ERROR = 40};
   State state = START;
   State prevState = START;

   Variable::Type variable;
   Options dOptions;
   Options cOptions;
   int index = 2;

   if(argv.size() == 2) {
      return false;
   }

   std::string downscaler = "nearestNeighbour";
   std::string calibrator = "";
   std::vector<Calibrator*> calibrators;
   while(true) {
      // std::cout << state << std::endl;
      if(state != ERROR)
         prevState = state;
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
         if(argv.size() <= index) {
            // -v but nothing after it
            state = ERROR;
         }
         else {
            variable = Variable::getType(argv[index]);
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else {
               state = ERROR;
            }
         }
      }
      else if(state == NEWVAR) {
         dOptions.addOption("variable", Variable::getTypeName(variable));
         Downscaler* d = Downscaler::getScheme(downscaler, variable, dOptions);
         VariableConfiguration varconf;
         varconf.variable = variable;
         varconf.downscaler = d;
         varconf.calibrators = calibrators;
         iSetup.variableConfigurations.push_back(varconf);

         // Reset to defaults
         downscaler = "nearestNeighbour";
         dOptions.clear();
         calibrators.clear();

         if(argv.size() <= index) {
            state = END;
         }
         else {
            state = VAR;
            index++;
         }
      }
      else if(state == DOWN) {
         if(argv.size() <= index) {
            // -d but nothing after it
            state = ERROR;
         }
         else {
            downscaler = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else {
               state = DOWNOPT;
            }
         }
      }
      else if(state == DOWNOPT) {
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else {
            // Process downscaler options
            if(argv.size() > index) {
               dOptions.addOption(argv[index]);
               index++;
            }
            /*
            if(argv.size() > index+1) {
               dOptions.addOption(argv[index], argv[index+1]);
               index++;
               index++;
            }
            */
            else {
               state = ERROR;
            }
         }
      }
      else if(state == CAL) {
         if(argv.size() <= index) {
            // -c but nothing after it
            state = ERROR;
         }
         else {
            calibrator = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWCAL;
            }
            else if(argv[index] == "-v") {
               state = NEWCAL;
            }
            else if(argv[index] == "-c") {
               state = NEWCAL;
            }
            else if(argv[index] == "-d") {
               state = NEWCAL;
            }
            else {
               state = CALOPT;
            }
         }
      }
      else if(state == CALOPT) {
         if(argv.size() <= index) {
            state = NEWCAL;
         }
         else if(argv[index] ==  "-c") {
            state = NEWCAL;
         }
         else if(argv[index] == "-v") {
            state = NEWCAL;
         }
         else {
            // Process calibrator options
            if(argv.size() > index) {
               cOptions.addOption(argv[index]);
               index++;
            }
            else {
               state = ERROR;
            }
         }
      }
      else if(state == NEWCAL) {
         cOptions.addOption("variable", Variable::getTypeName(variable));
         Calibrator* c = Calibrator::getScheme(calibrator, cOptions);
         calibrators.push_back(c);

         // Reset
         calibrator = "";
         cOptions.clear();
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] == "-d") {
            state = DOWN;
            index++;
         }
         else {
            state = ERROR;
         }
      }
      else if(state == END) {
         break;
      }
      else if(state == ERROR) {
         std::stringstream ss;
         ss << "Finite state machine entered error state. Previous state: " << prevState;
         Util::error(ss.str());
      }
   }
   return true;
}
