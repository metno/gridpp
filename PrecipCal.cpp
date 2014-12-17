#include <iostream>
#include <string>
#include "DataFile.h"
#include "ParameterFile.h"
#include "Calibration.h"
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
   std::string dataFile   = argv[1];
   std::string parameterFile = argv[2];
   std::string outputFile = argv[3];

   DataFile input(dataFile);
   DataFile output(outputFile);
   ParameterFile parameters(parameterFile);

   Calibration cal(parameters);
   cal.calibrate(input, output);

   output.write();
   return 0;
}
