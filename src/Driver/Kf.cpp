#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
#include "../Options.h"
#include "../SetupKf.h"
#include "../KalmanFilter.h"
void writeUsage() {
   std::cout << "Kalman filter program to create a bias file to be used by kriging in gridpp" << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  kalmanFilter [kfopt] -v var -o file -f file [-p output] [-dbin file] [-dbout file]" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << "   kfopt         Kalman filter options, see below." << std::endl;
   std::cout << "   -v var        One of the variables below." << std::endl;
   // std::cout << "   -d downscaler One of the downscalers below." << std::endl;
   std::cout << "   -o file       Observation file." << std::endl;
   std::cout << "   -f file       Forecast file." << std::endl;
   std::cout << "   -p output     Put biases into this file. One of the parameter formats below." << std::endl;
   std::cout << "   -dbin file    Input KF database." << std::endl;
   std::cout << "   -dbout file   Output KF database." << std::endl;
   std::cout << std::endl;

   std::cout << "Example:" << std::endl;
   std::cout << "   ./gridpp_kf time=0 -o testing/files/obs.txt type=text spatial=1\\" << std::endl;
   std::cout << "                      -f testing/files/10x10.nc\\" << std::endl;
   std::cout << "                      -v T -p text spatial=1 file=bias.txt " << std::endl;
   std::cout << std::endl;

   std::cout << "KF options:" << std::endl;
   std::cout << KalmanFilter::description();
   std::cout << std::endl;

   std::cout << "Obs/fcst files:" << std::endl;
   std::cout << File::getDescriptions();
   std::cout << std::endl;
   // std::cout << "Downscalers:" << std::endl;
   // std::cout << Downscaler::getDescriptions();
   // std::cout << std::endl;

   std::cout << "KF databases and output bias file:" << std::endl;
   std::cout << ParameterFile::getDescriptions();
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

   // Retrieve setup
   std::vector<std::string> args;
   for(int i = 1; i < argc; i++) {
      args.push_back(std::string(argv[i]));
   }
   SetupKf setup(args);

   int time = 0;
   setup.options.getValue("time", time);

   KalmanFilter kf(setup.variable, setup.options);
   kf.writeBiasFile(*setup.fcstFile, *setup.obsFile, time, setup.dbin, setup.dbout, setup.output);

}
