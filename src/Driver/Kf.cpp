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
#include "../KalmanFilter.h"
void writeUsage() {
   std::cout << "Kalman filter program to create a bias file to be used by kriging in gridpp" << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  kalmanFilter options" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << Util::formatDescription("data=required","Input text file with obs/forecast data. Must be in the following format:") << std::endl;
   std::cout << Util::formatDescription("","<time> <lat> <lon> <elev> <obs> <fcst>") << std::endl;
   std::cout << Util::formatDescription("current=undef","Current Kalman coefficients. If unspecified, initialize. Must have the format:") << std::endl;
   std::cout << Util::formatDescription("","<time> <lat> <lon> <elev> <pVarV> <kalmanGainVar> <varV> <p> <kalmanGain> <error> <lastError> <biasEstimate>") << std::endl;
   std::cout << Util::formatDescription("new=undef","New coeffcients stored in this file.") << std::endl;
   std::cout << Util::formatDescription("output=undef","Bias file, suitable for kriging.") << std::endl;
   std::cout << std::endl;
   std::cout << "Notes:" << std::endl;
   std::cout << "   - If 'input' contains stations that don't exist in currCoeffs, new parameters" << std::endl;
   std::cout << "     initialized." << std::endl;

   std::cout << std::endl;
}
int main(int argc, const char *argv[]) {
   if(argc <= 1) {
      writeUsage();
      return 0;
   }

   std::string dataFilename = "";
   std::string currentFilename = "";
   std::string newFilename = "";
   std::string outputFilename = "";
   Options options;
   for(int i = 1; i < argc; i++) {
      options.addOptions(std::string(argv[i]));
   }
   if(!options.getValue("data", dataFilename)) {
      Util::error("kalmanFilter: 'data' required.");
   }
   options.getValue("current", currentFilename);
   options.getValue("new", newFilename);
   options.getValue("output", outputFilename);

   KalmanFilter kf(Variable::T, 0.1);
   ParameterFileText dataFile(Options("file=" + dataFilename + " spatial=1"));
   ParameterFileText currFile(Options("file=" + currentFilename + " spatial=1"));
   ParameterFileText newFile(Options("file=" + newFilename + " spatial=1"));
   ParameterFileText outputFile(Options("file=" + outputFilename + " spatial=1"));

   std::vector<int>      times     = Util::combine(dataFile.getTimes(), currFile.getTimes());
   std::vector<Location> locations = Util::combine(dataFile.getLocations(), currFile.getLocations());
   for(int t = 0; t < times.size(); t++) {
      int time = times[t];
      for(int i = 0; i < locations.size(); i++) {
         const Location& location = locations[i];

         Parameters currData = dataFile.getParameters(time, location);
         // KF must still be run when obs/fcst are missing because it needs to update certain
         // coefficients.
         float obs      = Util::MV;
         float forecast = Util::MV;
         if(currData.size() == 2) {
            obs      = currData[0];
            forecast = currData[1];
         }

         std::vector<float> values;
         Parameters currKF = currFile.getParameters(time, location);
         Parameters newKF  = kf.update(obs, forecast, currKF);
         newFile.setParameters(newKF, time, location);
         if(outputFilename != "") {
            Parameters bias;
            std::vector<float> biasVector(1, newKF[7]);
            outputFile.setParameters(biasVector, time, location);
         }
      }
   }

   newFile.write();
   if(outputFilename != "")
      outputFile.write();
}
