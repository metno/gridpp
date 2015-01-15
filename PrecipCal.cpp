#include <iostream>
#include <string>
#include <string.h>
#include "File.h"
#include "ParameterFile.h"
#include "Calibrator/Calibrator.h"
#include "Util.h"
int main(int argc, const char *argv[]) {

   // Parse command line attributes
   if(argc < 3) {
      std::cout << "Calibrator of ensemble precipitation forecasts" << std::endl;
      std::cout << std::endl;
      std::cout << "usage:  precipCalibrator.exe file parameters [-v]" << std::endl;
      std::cout << std::endl;
      std::cout << "Arguments:" << std::endl;
      std::cout << "   file          Netcdf file with ECMWF ensemble data. Data in this file is overwritten." << std::endl;
      std::cout << "   parameters    Text file with parameters." << std::endl;
      std::cout << "   -v            Verbose. Show status messages." << std::endl;
      return 1;
   }
   double tStart = Util::clock();
   std::string dataFile      = argv[1];
   std::string parameterFile = argv[2];

   for(int i = 3; i < argc; i++) {
      if(strcmp(argv[i],"--debug"))
         Util::setShowStatus(true);
   }

   Util::setShowError(true);
   Util::setShowWarning(true);

   File file(dataFile);
   ParameterFile parameters(parameterFile);

   std::vector<Variable::Type> writableVariables;
   // Temperature calibration
   // CalibratorTemperature calTemp(parameters);
   // calTemp.calibrate(output);

   // Precip calibration
   CalibratorZaga calPrecip(parameters, Variable::Precip);
   calPrecip.setFracThreshold(0.5);
   calPrecip.calibrate(file);

   // Cloud calibration
   CalibratorCloud calCloud(parameters);
   calCloud.calibrate(file);
   writableVariables.push_back(Variable::Cloud);

   CalibratorAccumulate calPrecipAcc(parameters, Variable::Precip, Variable::PrecipAcc);
   calPrecipAcc.calibrate(file);
   writableVariables.push_back(Variable::PrecipAcc);

   // Wind calibration
   // CalibratorWind calWind(parameters);
   // calWind.calibrate(output);

   // Output
   double t3 = Util::clock();
   file.write(writableVariables);
   double tEnd = Util::clock();
   std::cout << "Calibrator time: " << t3 - tStart << std::endl;
   std::cout << "Writing time: " << tEnd - t3 << std::endl;
   std::cout << "Total time: " << tEnd - tStart << std::endl;
   return 0;
}
