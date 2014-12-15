#include <iostream>
#include <string>
#include "DataFile.h"
#include "ParameterFile.h"
#include "Calibration.h"
#include <ctime>
int main(int argc, const char *argv[]) {

   // Parse command line attributes
   if(argc != 4) {
      std::cout << "Usage:" << std::endl;
      std::cout << "   precipCalibration.exe <dataFile> <parameterFile> <outputFile>" << std::endl;
      std::cout << "   - dataFile: Netcdf file with Arome data" << std::endl;
      std::cout << "   - parameterFile: Text file with training coefficients" << std::endl;
      std::cout << "   - outputFile: Name of file to write factors" << std::endl;
      return 1;
   }
   std::string dataFile   = argv[1];
   std::string parameterFile = argv[2];
   std::string outputFile = argv[3];

   // Get list of coefficients, one set for each site
   ParameterFile parameters(parameterFile);

   // Set up data files
   DataFile input(dataFile);
   DataFile output(outputFile);

   Calibration cal(parameters);
   cal.calibrate(input, output);

   output.write();
   return 0;
}
