#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
#include "../Options.h"
#include "../SetupTrain.h"
#include "../TrainingData.h"
void writeUsage() {
   std::cout << "Trains a calibration method in gridpp" << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  gridpp_train dataFile -p parameters [options] -c calibrator [options] -v variable" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << Util::formatDescription("dataFile=required","Input text file with obs/forecast data. Must be in the following format:") << std::endl;
   std::cout << Util::formatDescription("","<time> <lat> <lon> <elev> <obs> <ens1> <ens2> ... <ensN>") << std::endl;
   std::cout << Util::formatDescription("-p parameters","Put parameters into this output format.") << std::endl;
   std::cout << Util::formatDescription("-c calibrator","Train this calibration method") << std::endl;
   std::cout << Util::formatDescription("-v variable","Train this variable") << std::endl;
   std::cout << std::endl;
   std::cout << "Example:" << std::endl;
   std::cout << "   ./gridpp_train testing/files/training.txt -p text file=parameters.txt -c regression -v T" << std::endl;
   std::cout << std::endl;
   std::cout << "Variables:" << std::endl;
   std::cout << Variable::getDescriptions();
   std::cout << std::endl;
   std::cout << "Calibrators with options (and default values):" << std::endl;
   std::cout << Calibrator::getDescriptions();
   std::cout << std::endl;
   std::cout << "Parameter formats with Options (and default values)" << std::endl;
   std::cout << ParameterFile::getDescriptions();
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
   SetupTrain setup(args);

   File* forecast = setup.forecasts[0];
   File* observation = setup.observations[0];
   int nLat = observation->getNumLat();
   int nLon = observation->getNumLon();
   int nEns = observation->getNumEns();
   int nTime = observation->getNumTime();
   vec2 lats = observation->getLats();
   vec2 lons = observation->getLons();
   vec2 elevs = observation->getElevs();

   Variable::Type variable = Variable::T;
   std::vector<double> offsets = forecast->getTimes();
   int D = setup.forecasts.size();

   std::cout << "Number of files: " << D << std::endl;

   for(int t = 0; t < offsets.size(); t++) {
      // Loop over all files
      std::vector<FieldPtr> ffields;
      std::vector<FieldPtr> ofields;
      for(int d = 0; d < D; d++) {
         ffields.push_back(setup.forecasts[d]->getField(variable, t));
         ofields.push_back(setup.observations[d]->getField(variable, t));
      }

      vec2Int I, J;
      Downscaler::getNearestNeighbour(*forecast, *observation, I, J);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int offset = t;
            double timeStart = Util::clock();

            // Arrange data
            std::vector<ObsEns> data;
            for(int d = 0; d < D; d++){
               // Downscaling (currently nearest neighbour)
               float obs = (*ofields[d])(I[i][j],J[i][j],0);
               const Ens& ens = (*ffields[d])(i,j);
               ObsEns obsens(obs, ens);
               data.push_back(obsens);
            }

            Parameters par = setup.method->train(data);
            double timeEnd = Util::clock();
            std::cout << "Time: " << timeEnd - timeStart << std::endl;
            std::cout << "Calculating parameters for offset=" << i << ": ";
            for(int k = 0; k < par.size(); k++) {
               std::cout << " " << par[k];
            }
            std::cout << std::endl;
            setup.output->setParameters(par, offset, Location(Util::MV, Util::MV, Util::MV));
         }
      }
   }

   setup.output->write();
}
