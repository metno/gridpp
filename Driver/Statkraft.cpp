#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
int main(int argc, const char *argv[]) {

   // Parse command line attributes
   if(argc < 3) {
      std::cout << "Calibrator of ensemble precipitation forecasts" << std::endl;
      std::cout << std::endl;
      std::cout << "usage:  statkraft.exe input output parameters [-v]" << std::endl;
      std::cout << std::endl;
      std::cout << "Arguments:" << std::endl;
      std::cout << "   input         Netcdf file with AROME data." << std::endl;
      std::cout << "   output        Netcdf file with AROME data." << std::endl;
      std::cout << "   parameters    Text file with parameters." << std::endl;
      std::cout << "   -v            Verbose. Show status messages." << std::endl;
      return 1;
   }
   double tStart = Util::clock();
   std::string inputFile     = argv[1];
   std::string outputFile    = argv[2];
   // std::string parameterFile = argv[3];

   for(int i = 3; i < argc; i++) {
      if(strcmp(argv[i],"--debug"))
         Util::setShowStatus(true);
   }

   Util::setShowError(true);
   Util::setShowWarning(true);

   const FileArome ifile(inputFile);
   FileArome ofile(outputFile);
   // ParameterFileRegion parameters(parameterFile);

   std::vector<Variable::Type> writableVariables;

   // Downscaling
   double t3 = Util::clock();
   DownscalerGradient downscaler(Variable::T);
   downscaler.downscale(ifile, ofile);
   writableVariables.push_back(Variable::T);

   // Wind calibration
   // CalibratorWind calWind(parameters);
   // calWind.calibrate(ofile);

   // Output
   double t4 = Util::clock();
   ofile.write(writableVariables);
   double tEnd = Util::clock();
   std::cout << "Downscaler time: " << t3 - tStart << std::endl;
   std::cout << "Calibrator time: " << t4 - tStart << std::endl;
   std::cout << "Writing time: " << tEnd - t4 << std::endl;
   std::cout << "Total time: " << tEnd - tStart << std::endl;
   return 0;
}
