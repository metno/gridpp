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
   std::cout << "Inputs/Outputs:" << std::endl;
   std::cout << KalmanFilter::description();
   std::cout << std::endl;
}
int main(int argc, const char *argv[]) {
   if(argc <= 1) {
      writeUsage();
      return 0;
   }

   std::string fcstFilename = "";
   std::string obsFilename = "";
   std::string dbFilename = "";
   std::string outputFilename = "";
   Options options;
   for(int i = 1; i < argc; i++) {
      options.addOptions(std::string(argv[i]));
   }
   if(!options.getValue("fcst", fcstFilename)) {
      Util::error("kalmanFilter: 'fcst' required.");
   }
   options.getValue("obs", obsFilename);
   int startTime = 0;
   int endTime = 0;
   options.getValue("startTime", startTime);
   options.getValue("endTime", endTime);

   Options kfOptions = options;

   KalmanFilter kf(Variable::T, kfOptions);
   FilePoint fcstFile = FilePoint(fcstFilename, Options("lat=5 lon=5 elev=1000 ens=1"));

   // Downscale to observations

   FilePoint obsFile = FilePoint(obsFilename, Options("lat=5 lon=5 elev=1000 ens=1"));
   ParameterFileText* dbIn = NULL;
   if(options.getValue("dbin", dbFilename)) {
      dbIn = new ParameterFileText(Options("file=" + dbFilename + " spatial=1"));
   }
   ParameterFileText* dbOut = NULL;
   if(options.getValue("dbout", dbFilename)) {
      dbOut = new ParameterFileText(Options("file=" + dbFilename + " spatial=1"));
   }
   ParameterFileText* outputFile = NULL;
   if(options.getValue("output", outputFilename)) {
      outputFile = new ParameterFileText(Options("file=" + outputFilename + " spatial=1"));
   }

   kf.writeBiasFile(fcstFile, obsFile, 20150101, startTime, endTime, dbIn, dbOut, outputFile);

   if(dbIn != NULL)
      delete dbIn;
   if(dbOut != NULL)
      delete dbOut;
   if(outputFile != NULL)
      delete outputFile;
}
