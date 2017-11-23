#include "Setup.h"
#include "File/File.h"
#include "Calibrator/Calibrator.h"
#include "Downscaler/Downscaler.h"

Setup::Setup(const std::vector<std::string>& argv) {

   std::string inputFilename = "";
   std::string outputFilename = "";

   // Process input/output filenames and options
   int index = 0;
   while(index < argv.size()) {
      std::string arg = argv[index];
      if(inputFilename == "") {
         inputFilename = arg;
      }
      else if(outputFilename == "") {
         if(Util::hasChar(arg, '=')) {
            inputOptions.addOptions(arg);
         }
         else {
            outputFilename = arg;
         }
      }
      else {
         if(Util::hasChar(arg, '=')) {
            outputOptions.addOptions(arg);
         }
         else {
            break;
         }
      }
      index++;
   }
   std::vector<std::string> inputFilenames = Util::glob(inputFilename);
   std::vector<std::string> outputFilenames = Util::glob(outputFilename);
   if(inputFilenames.size() != outputFilenames.size()) {
      std::stringstream ss;
      ss << "Unequal number of input (" << inputFilenames.size() << ") and output ("
         << outputFilenames.size() << ")";
      Util::error(ss.str());
   }
   if(inputFilenames.size() == 0) {
      Util::error("No valid input files");
   }
   if(outputFilenames.size() == 0) {
      Util::error("No valid output files");
   }

   for(int f = 0; f < outputFilenames.size(); f++) {

      File* outputFile;
      if(!hasFile(outputFilenames[f])) {
         outputFile  = File::getScheme(outputFilenames[f], outputOptions, false);
         if(outputFile == NULL) {
            Util::error("File '" + outputFilenames[f] + " is invalid");
         }
         mFileMap[outputFilenames[f]] = outputFile;
      }
      else {
         outputFile = mFileMap[outputFilenames[f]];
      }
      outputFiles.push_back(outputFile);

      File* inputFile;
      if(!hasFile(inputFilenames[f])) {
         inputFile  = File::getScheme(inputFilenames[f], inputOptions, true);
         if(inputFile == NULL) {
            Util::error("File '" + inputFilenames[f] + " is invalid");
         }
         mFileMap[inputFilenames[f]] = inputFile;
      }
      else {
         inputFile = mFileMap[inputFilenames[f]];
      }
      inputFiles.push_back(inputFile);
   }

   // Implement a finite state machine
   enum State {START = 0, VAR = 1, VAROPT = 2, NEWVAR = 3, DOWN = 10, DOWNOPT = 15, CAL = 20, NEWCAL = 22, CALOPT = 25, PARDOWN = 30, PAROPTDOWN = 35, PARCAL = 40, PAROPTCAL = 45, END = 90, ERROR = 100};
   State state = START;
   State prevState = START;

   std::string variableName = "";
   std::string variableGridppName = "";
   Options vOptions;
   Options dOptions;
   Options cOptions;
   Options pOptions;
   std::string errorMessage = "";

   std::string downscaler = defaultDownscaler();
   std::string calibrator = "";
   std::string parameterFile = "";
   std::vector<std::string> calibratorNames;
   std::vector<ParameterFile*> parameterFileCalibrators;
   ParameterFile* parameterFileDownscaler = NULL;
   while(true) {
      // std::cout << state << std::endl;
      if(state != ERROR)
         prevState = state;
      if(state == START) {
         if(index < argv.size() && argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else {
            errorMessage = "No variables defined";
            state = ERROR;
         }
      }
      else if(state == VAR) {
         if(argv.size() <= index) {
            // -v but nothing after it
            errorMessage = "No variable after '-v'";
            state = ERROR;
         }
         else {
            variableGridppName = argv[index];
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
            else if(argv[index] == "-p") {
               errorMessage = "-p must be after a -d or -c";
               state = ERROR;
            }
            else {
               state = VAROPT;
            }
         }
      }
      else if(state == VAROPT) {
         if(argv.size() <= index) {
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
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] == "-p") {
            errorMessage = "-p must be after a -d or -c";
            state = ERROR;
         }
         else {
            // Process variable options
            vOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == NEWVAR) {
         // Check that we haven't added the variable before
         bool alreadyExists = false;
         for(int i = 0; i < variableConfigurations.size(); i++) {
            if(variableConfigurations[i].inputVariable.getName() == variableGridppName)
               alreadyExists = true;
         }

         if(!alreadyExists) {
            VariableConfiguration varconf;
            varconf.parameterFileDownscaler = parameterFileDownscaler;
            // TODO: 0
            bool found = inputFiles[0]->getVariable(Variable::getType(variableGridppName), varconf.inputVariable);
            if(!found) {
               std::stringstream ss;
               ss << "Input format does not define variable of type '" << variableGridppName << "'";
               Util::error(ss.str());
            }
            found = outputFiles[0]->getVariable(Variable::getType(variableGridppName), varconf.outputVariable);
            if(!found) {
               std::stringstream ss;
               ss << "Output format does not define variable of type '" << variableGridppName << "'";
               Util::error(ss.str());
            }
            std::vector<Calibrator*> calibrators;
            for(int c = 0; c < calibratorNames.size(); c++) {
               Calibrator* calibrator = Calibrator::getScheme(calibratorNames[c], varconf.outputVariable, cOptions);
               calibrators.push_back(calibrator);
            }
            varconf.calibrators = calibrators;
            varconf.parameterFileCalibrators = parameterFileCalibrators;
            varconf.variableOptions = vOptions;
            Downscaler* d = Downscaler::getScheme(downscaler, varconf.inputVariable, varconf.outputVariable, dOptions);
            varconf.downscaler = d;
            variableConfigurations.push_back(varconf);
         }
         else {
            std:: stringstream ss;
            ss << "Variable '" << variableGridppName << "' already read. Using first instance.";
            Util::warning(ss.str());
         }

         // Reset to defaults
         vOptions.clear();
         downscaler = defaultDownscaler();
         parameterFileDownscaler = NULL;
         dOptions.clear();
         calibratorNames.clear();
         parameterFileCalibrators.clear();

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
            errorMessage = "No downscaler after '-d'";
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
            else if(argv[index] == "-d") {
               // Two downscalers defined for one variable
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-p") {
               state = PARDOWN;
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
         else if(argv[index] == "-p") {
            state = PARDOWN;
         }
         else {
            // Process downscaler options
            dOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == PARDOWN) {
         if(argv.size() <= index) {
            // -p but nothing after it
            errorMessage = "No parameter file after '-p'";
            state = ERROR;
         }
         else {
            parameterFile = argv[index];
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
            else if(argv[index] == "-d") {
               // Two downscalers defined for one variable
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-p") {
               // Two parameter files defined for one variable
               errorMessage = "Two or more -p used for one downscaler";
               state = ERROR;
               index++;
            }
            else {
               state = PAROPTDOWN;
            }
         }
      }
      else if(state == PAROPTDOWN) {
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
         else if(argv[index] == "-p") {
            errorMessage = "Two or more -p used for one calibrator";
            state = ERROR;
         }
         else {
            // Process parameter file options
            pOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == CAL) {
         if(argv.size() <= index) {
            // -c but nothing after it
            state = ERROR;
            errorMessage = "No calibrator after '-c'";
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
            else if(argv[index] == "-p") {
               state = PARCAL;
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
         else if(argv[index] == "-d") {
            state = NEWCAL;
         }
         else if(argv[index] == "-p") {
            state = PARCAL;
         }
         else {
            // Process calibrator options
            cOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == PARCAL) {
         if(argv.size() <= index) {
            // -p but nothing after it
            state = ERROR;
            errorMessage = "No parameter file after '-p'";
         }
         else {
            parameterFile = argv[index];
            if(argv.size() <= index) {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-c") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-d") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-p") {
               state = PARCAL;
               index++;
            }
            else {
               state = PAROPTCAL;
            }
         }
      }
      else if(state == PAROPTCAL) {
         if(argv.size() <= index) {
            state = NEWCAL;
         }
         else if(argv[index] ==  "-c") {
            state = NEWCAL;
         }
         else if(argv[index] == "-v") {
            state = NEWCAL;
         }
         else if(argv[index] == "-d") {
            state = NEWCAL;
         }
         else if(argv[index] == "-p") {
            state = PARCAL;
         }
         else {
            // Process calibrator options
            pOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == NEWCAL) {
         // We do not need to check that the same calibrator has been added for this variable
         // since this is perfectly fine (e.g. smoothing twice).

         ParameterFile* p = NULL;
         if(parameterFile != "") {
            std::string schemeName;
            if(!pOptions.getValue("type", schemeName)) {
               state = ERROR;
               errorMessage = "Parameter file missing 'type': " + pOptions.toString();
            }
            else {
               pOptions.addOption("file", parameterFile);
               p = ParameterFile::getScheme(schemeName, pOptions);
               if(!p->isReadable()) {
                  state = ERROR;
                  errorMessage = "Could not open parameter file: " + pOptions.toString();
               }
            }
         }
         if(state != ERROR) {
            calibratorNames.push_back(calibrator);
            parameterFileCalibrators.push_back(p);

            // Reset
            calibrator = "";
            parameterFile = "";
            cOptions.clear();
            pOptions.clear();
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
               errorMessage = "No recognized option after '-c calibrator'";
            }
         }
      }
      else if(state == END) {
         break;
      }
      else if(state == ERROR) {
         std::stringstream ss;
         ss << "Invalid command line arguments: " << errorMessage << ".";
         Util::error(ss.str());
      }
   }
}
Setup::~Setup() {
   std::map<std::string, File*>::const_iterator it = mFileMap.begin();
   for(it = mFileMap.begin(); it != mFileMap.end(); it++)
      delete it->second;

   for(int i = 0; i < variableConfigurations.size(); i++) {
      delete variableConfigurations[i].downscaler;
      for(int c = 0; c < variableConfigurations[i].calibrators.size(); c++) {
         delete variableConfigurations[i].calibrators[c];
      }
   }
}
std::string Setup::defaultDownscaler() {
   return "nearestNeighbour";
}

bool Setup::hasFile(std::string iFilename) const {
   std::map<std::string, File*>::const_iterator it = mFileMap.find(iFilename);
   return it != mFileMap.end();
}
