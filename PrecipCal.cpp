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
      std::cout << "   output        Netcdf file to write output to. Must already exist and have \n"
                << "                 the same dimensions as 'input'. 'output' may be the same file as 'input'." << std::endl;
      return 1;
   }
   double tStart = Util::clock();
   std::string dataFile      = argv[1];
   std::string parameterFile = argv[2];
   std::string outputFile    = argv[3];

   Util::setShowError(true);
   Util::setShowWarning(true);
   Util::setShowStatus(true);

   DataFile input(dataFile);
   DataFile output(outputFile);
   ParameterFile parameters(parameterFile);

   // Checks
   if(!input.hasSameDimensions(output)) {
      std::stringstream ss;
      ss << "'" << input.getFilename() << "' " << input.getDimenionString() << " does not have the same dimensions"
         << " as '" << output.getFilename() << "' " << output.getDimenionString();
      Util::error(ss.str());
   }

   // Calibration
   CalibrationPrecip cal(parameters);
   cal.setFracThreshold(0.5);
   cal.calibrate(input, output);

   // Output
   double t3 = Util::clock();
   output.write();
   double tEnd = Util::clock();
   std::cout << "Calibration time: " << t3 - tStart << std::endl;
   std::cout << "Writing time: " << tEnd - t3 << std::endl;
   std::cout << "Total time: " << tEnd - tStart << std::endl;
   return 0;
}
