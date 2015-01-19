#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
float getGradient(int iDate);
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
   std::string WparameterFilename = "test.txt";//argv[3];

   for(int i = 3; i < argc; i++) {
      if(strcmp(argv[i],"--debug"))
         Util::setShowStatus(true);
   }

   Util::setShowError(true);
   Util::setShowWarning(true);

   const FileArome ifile(inputFile, true);
   int date = ifile.getDate();
   FileArome ofile(outputFile);

   std::vector<Variable::Type> writableVariables;

   /////////////////
   // Temperature //
   /////////////////
#if 1
   // Downscaling
   double tt0 = Util::clock();
   float gradient = getGradient(date);
   std::cout << "Temperature gradient: " << gradient << " " << date << std::endl;
   DownscalerGradient Tdownscaler(Variable::T);
   Tdownscaler.setConstantGradient(gradient);
   Tdownscaler.downscale(ifile, ofile);
   writableVariables.push_back(Variable::T);
   double tt1 = Util::clock();
   std::cout << "Temperature time: " << tt1 - tt0 << std::endl;

   //////////
   // Wind //
   //////////
   // TODO: Can you use the gradient approach on U and V separately?
   double tw0 = Util::clock();
   ofile.initNewVariable(Variable::U);
   ofile.initNewVariable(Variable::V);
   DownscalerGradient Udownscaler(Variable::U);
   Udownscaler.setSearchRadius(7);
   Udownscaler.downscale(ifile, ofile);
   DownscalerGradient Vdownscaler(Variable::V);
   Vdownscaler.setSearchRadius(7);
   Vdownscaler.downscale(ifile, ofile);

   if(0) {
      ParameterFileRegion Wparameters(WparameterFilename);
      CalibratorWind Wcalibrator(Wparameters);
      Wcalibrator.calibrate(ofile);
   }
   writableVariables.push_back(Variable::U);
   writableVariables.push_back(Variable::V);
   double tw1 = Util::clock();
   std::cout << "Wind time:        " << tw1 - tw0 << std::endl;
#endif

   ////////////
   // Precip //
   ////////////
   double tp0 = Util::clock();
   ofile.initNewVariable(Variable::PrecipAcc);
   DownscalerSmart Pdownscaler(Variable::Precip);
   Pdownscaler.setSearchRadius(5);
   Pdownscaler.setNumSmart(5);
   Pdownscaler.downscale(ifile, ofile);

   CalibratorAccumulate Paccumulate(Variable::Precip, Variable::PrecipAcc);
   Paccumulate.calibrate(ofile);
   writableVariables.push_back(Variable::Precip);
   writableVariables.push_back(Variable::PrecipAcc);
   double tp1 = Util::clock();
   std::cout << "Precip time:      " << tp1 - tp0 << std::endl;

   // Output
   double t4 = Util::clock();
   ofile.write(writableVariables);
   double tEnd = Util::clock();
   std::cout << "Writing time:     " << tEnd - t4 << std::endl;
   std::cout << "Total time:       " << tEnd - tStart << std::endl;
   return 0;
}

float getGradient(int iDate) {
   int month = (iDate % 10000) / 100;
   if(month <= 3 || month > 9) {
      // Winter
      return -0.007;
   }
   else {
      // Summer
      return -0.006;
   }
}
