#include <iostream>
#include <string>
#include "DataFile.h"
#include "ParameterFile.h"
#include "Calibration.h"
#include "Util.h"
int main(int argc, const char *argv[]) {

   // Parse command line attributes
   if(argc != 4) {
      std::cout << "Calibration of ensemble precipitation forecasts" << std::endl;
      std::cout << std::endl;
      std::cout << "usage:  precipCalibration.exe input parameters output" << std::endl;
      std::cout << std::endl;
      std::cout << "Arguments:" << std::endl;
      std::cout << "   input         Netcdf file with ECMWF ensemble data" << std::endl;
      std::cout << "   parameters    Text file with parameters" << std::endl;
      std::cout << "   output        Netcdf file to write output to" << std::endl;
      return 1;
   }
   double tStart = Util::clock();
   std::string dataFile   = argv[1];
   std::string parameterFile = argv[2];
   std::string outputFile = argv[3];

   DataFile input(dataFile);
   DataFile output(outputFile);
   ParameterFile parameters(parameterFile);

   CalibrationPrecip cal(parameters);
   cal.calibrate(input, output);

   double t3 = Util::clock();
   output.write();
   double tEnd = Util::clock();
   std::cout << "Calibration time: " << t3 - tStart << std::endl;
   std::cout << "Writing time: " << tEnd - t3 << std::endl;
   std::cout << "Total time: " << tEnd - tStart << std::endl;
   return 0;
}
