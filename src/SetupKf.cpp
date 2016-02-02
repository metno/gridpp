#include "SetupKf.h"
#include "File/File.h"
#include "Calibrator/Calibrator.h"
#include "Downscaler/Downscaler.h"
#include "TrainingData.h"

SetupKf::SetupKf(const std::vector<std::string>& argv) {

   if(argv.size() <= 0) {
      Util::error("No inputs");
   }

   // Set initial non-working values
   variable = Variable::None;
   output = NULL;
   fcstFile = NULL;
   obsFile = NULL;
   dbin = NULL;
   dbout = NULL;
   downscaler = NULL; // Downscaler::getScheme("nearestNeighbour", Variable::T, Options());

   // Process obs/fcst filenames and options
   Options obsOptions, fcstOptions;
   std::string obsFilename, fcstFilename;
   int index = 0;
   while(index < argv.size()) {
      if(argv[index][0] == '-')
         break;
      options.addOptions(argv[index]);
      index++;
   }

   while(index < argv.size()) {
      if(argv[index] == "-v") {
         index++;
         variable = Variable::getType(argv[index]);
         index++;
      }
      else if(argv[index] == "-o") {
         index++;
         std::string name = argv[index];
         Options opt;
         index++;
         while(index < argv.size()) {
            if(argv[index][0] == '-') 
               break;
            opt.addOptions(argv[index]);
            index++;
         }
         obsFile = File::getScheme(name, opt, true);
      }
      else if(argv[index] == "-f") {
         index++;
         std::string name = argv[index];
         Options opt;
         index++;
         while(index < argv.size()) {
            if(argv[index][0] == '-') 
               break;
            opt.addOptions(argv[index]);
            index++;
         }
         fcstFile = File::getScheme(name, opt, true);
      }
      else if(argv[index] == "-d") {
         index++;
         std::string downscalerName = argv[index];
         Options dOptions;
         index++;
         while(index < argv.size()) {
            if(argv[index][0] == '-') 
               break;
            dOptions.addOptions(argv[index]);
            index++;
         }
         if(variable == Variable::None) {
            Util::error("-d before -v");
         }
         downscaler = Downscaler::getScheme(downscalerName, variable, dOptions);
      }
      else if(argv[index] == "-dbin") {
         index++;
         std::string name = argv[index];
         Options opt;
         index++;
         while(index < argv.size()) {
            if(argv[index][0] == '-') 
               break;
            opt.addOptions(argv[index]);
            index++;
         }
         dbin = ParameterFile::getScheme(name, opt);
      }
      else if(argv[index] == "-dbout") {
         index++;
         std::string name = argv[index];
         Options opt;
         index++;
         while(index < argv.size()) {
            if(argv[index][0] == '-') 
               break;
            opt.addOptions(argv[index]);
            index++;
         }
         dbout = ParameterFile::getScheme(name, opt);
      }
      else if(argv[index] == "-p") {
         index++;
         std::string name = argv[index];
         index++;
         Options opt;
         while(index < argv.size()) {
            if(argv[index][0] == '-') 
               break;
            opt.addOptions(argv[index]);
            index++;
         }
         output = ParameterFile::getScheme(name, opt, true);
      }
      else {
         Util::error("Could not recognize '" + argv[index] + "'");
      }
   }
   if(fcstFile == NULL)
      Util::error("Missing fcstFile");
   if(obsFile == NULL)
      Util::error("Missing obsFile");
   // if(downscaler == NULL)
   //    Util::error("Missing downscaler");
   if(variable == Variable::None)
      Util::error("Missing variable");
   if(output == NULL)
      Util::error("Missing output");
}
SetupKf::~SetupKf() {
   delete obsFile;
   delete fcstFile;
   delete dbin;
   delete dbout;
   delete output;
   delete downscaler;
}
