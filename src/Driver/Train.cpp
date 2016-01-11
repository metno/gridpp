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
   Util::setShowError(true);
   Util::setShowWarning(true);
   Util::setShowStatus(false);
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
   int nLat = forecast->getNumLat();
   int nLon = forecast->getNumLon();
   int nEns = forecast->getNumEns();
   int nTime = forecast->getNumTime();
   vec2 lats = forecast->getLats();
   vec2 lons = forecast->getLons();
   vec2 elevs = forecast->getElevs();

   Variable::Type variable = Variable::T;
   std::vector<double> offsets = forecast->getTimes();
   int D = setup.forecasts.size();

   std::cout << "Number of files: " << D << std::endl;
   std::cout << "Number of lats: " << nLat << std::endl;
   std::cout << "Number of lons: " << nLon << std::endl;
   std::cout << "Number of times: " << offsets.size() << std::endl;
   // assert(observation->getNumTime() <= forecast->getNumTime());

   // figure out which files match
   std::vector<int> indexF;
   std::vector<int> indexO;
   for(int f = 0; f < setup.forecasts.size(); f++) {
      int io = Util::MV;
      double timeF = setup.forecasts[f]->getReferenceTime();
      for(int o = 0; o < setup.observations.size(); o++) {
         double timeO = setup.observations[f]->getReferenceTime();
         if(timeF == timeO) {
            io = o;
         }
      }
      if(Util::isValid(io)) {
         indexF.push_back(f);
         indexO.push_back(io);
      }
   }

   int validD = indexF.size();
   std::cout << "Number of valid days: " << validD << std::endl;
   if(validD == 0) {
      return -1;
   }

   vec2Int I, J;
   Downscaler::getNearestNeighbour(*observation, *forecast, I, J);

   for(int t = 0; t < 10; t++) { // offsets.size(); t++) {
      std::cout << t << std::endl;
      // Loop over all files
      std::vector<FieldPtr> ffields;
      std::vector<FieldPtr> ofields;
      for(int d = 0; d < validD; d++) {
         ffields.push_back(setup.forecasts[indexF[d]]->getField(variable, t));
         ofields.push_back(setup.observations[indexO[d]]->getField(variable, t));
      }

      // #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int offset = t;
            double timeStart = Util::clock();
            float lat = lats[i][j];
            float lon = lons[i][j];
            float elev = elevs[i][j];

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
            // std::cout << "Time: " << timeEnd - timeStart << std::endl;
            // std::cout << "Calculating parameters for offset=" << i << ": ";
            // for(int k = 0; k < par.size(); k++) {
            //    std::cout << " " << par[k];
            // }
            // std::cout << std::endl;
            setup.output->setParameters(par, offset, Location(lat, lon, elev));
         }
      }
   }

   setup.output->write();
}
