#include "SetupTrain.h"
#include "File/File.h"
#include "Calibrator/Calibrator.h"
#include "Downscaler/Downscaler.h"
#include "TrainingData.h"

SetupTrain::SetupTrain(const std::vector<std::string>& argv) {

   if(argv.size() <= 0) {
      Util::error("Missing training data file");
   }

   Variable::Type variable = Variable::T;

   // Implement a finite state machine
   enum State {START = 0, OUTPUT = 1, OUTPUTOPT = 2, METHOD = 3, METHODOPT = 10, VAR = 20, END = 90, ERROR = 100};
   State state = START;
   State prevState = START;

   std::string dataFilename = argv[0];
   trainingData = new TrainingData(dataFilename);
   if(argv.size() == 1) {
      Util::error("No method defined");
   }

   Options oOptions;
   Options mOptions;
   std::string outputType = "";
   std::string methodType = "";
   int index = 1;
   std::string errorMessage = "";
   while(true) {
      // std::cout << state << std::endl;
      if(state != ERROR)
         prevState = state;
      if(state == START) {
         if(index < argv.size() && argv[index] == "-o") {
            state = OUTPUT;
            index++;
         }
         else if(index < argv.size() && argv[index] == "-c") {
            state = METHOD;
            index++;
         }
         else if(index < argv.size() && argv[index] == "-v") {
            state = METHOD;
            index++;
         }
         else {
            errorMessage = "'-o', '-v', or '-c' reqired";
            state = ERROR;
         }
      }
      else if(state == OUTPUT) {
         if(argv.size() <= index) {
            // -o but nothing after it
            errorMessage = "Nothing after '-o'";
            state = ERROR;
         }
         else {
            outputType = argv[index];
            index++;
            if(argv.size() <= index) {
               state = END;
            }
            else if(argv[index] == "-o") {
               state = ERROR;
               errorMessage = "Duplicate '-o'";
            }
            else if(argv[index] == "-c") {
               state = METHOD;
               index++;
            }
            else if(argv[index] == "-v") {
               state = VAR;
               index++;
            }
            else {
               state = OUTPUTOPT;
            }
         }
      }
      else if(state == OUTPUTOPT) {
         if(argv.size() <= index) {
            state = END;
         }
         else if(argv[index] == "-c") {
            state = METHOD;
            index++;
         }
         else if(argv[index] == "-o") {
            state = ERROR;
            errorMessage = "Duplicate '-o'";
         }
         else if(argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else {
            // Process output options
            oOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == METHOD) {
         if(argv.size() <= index) {
            // -c but nothing after it
            errorMessage = "Nothing after '-c'";
            state = ERROR;
         }
         else {
            methodType = argv[index];
            index++;
            if(argv.size() <= index) {
               state = END;
            }
            else if(argv[index] == "-c") {
               state = ERROR;
               errorMessage = "Duplicate '-c'";
            }
            else if(argv[index] == "-o") {
               state = OUTPUT;
               index++;
            }
            else if(argv[index] == "-v") {
               state = VAR;
               index++;
            }
            else {
               state = METHODOPT;
            }
         }
      }
      else if(state == METHODOPT) {
         if(argv.size() <= index) {
            state = END;
         }
         else if(argv[index] == "-o") {
            state = OUTPUT;
            index++;
         }
         else if(argv[index] == "-c") {
            state = ERROR;
            errorMessage = "Duplicate '-c'";
         }
         else if(argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else {
            // Process method options
            mOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == VAR) {
         if(argv.size() <= index) {
            // -v but nothing after it
            errorMessage = "Nothing after '-v'";
            state = ERROR;
         }
         else {
            variable = Variable::getType(argv[index]);
            index++;
            if(argv.size() <= index) {
               state = END;
            }
            else if(argv[index] == "-c") {
               state = METHOD;
            }
            else if(argv[index] == "-o") {
               state = OUTPUT;
            }
            else if(argv[index] == "-v") {
               state = ERROR;
               errorMessage = "Duplicate '-v'";
            }
         }
      }
      else if(state == END) {
         mOptions.addOption("variable", Variable::getTypeName(variable));
         method = Calibrator::getScheme(methodType, mOptions);
         output = ParameterFile::getScheme(outputType, oOptions);
         break;
      }
      else if(state == ERROR) {
         std::stringstream ss;
         ss << "Could not understand command line arguments: " << errorMessage << ".";
         Util::error(ss.str());
      }
      else {
         abort();
      }
   }
   if(method == NULL) {
      Util::error("Missing method");
   }
   if(output == NULL) {
      Util::error("Missing output");
   }
}
SetupTrain::~SetupTrain() {
   delete trainingData;
   delete output;
   delete method;
}
