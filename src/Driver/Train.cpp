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
   std::cout << "Trains a calibration method in gridpp using observations and forecast files. The resulting parameter file is on the same grid as the forecasts and observations will be interpolated to the forecast grid using nearest neighbour." << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  gridpp_train obsFiles [options] fcstFiles [options] -p parameters [options] -c calibrator [options] -v variable" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << Util::formatDescription("obsFiles=required","Input file(s) with observations. Can be a regular expression if the whole string is surrounded by \" \". ") << std::endl;
   std::cout << Util::formatDescription("fcstFiles=required","Input file(s) with forecasts. Can be a regular expression if the whole string is surrounded by \" \". ") << std::endl;
   std::cout << Util::formatDescription("-p parameters","Put parameters into this output format.") << std::endl;
   std::cout << Util::formatDescription("-c calibrator","Train this calibration method") << std::endl;
   std::cout << Util::formatDescription("-v variable","Train this variable") << std::endl;
   std::cout << std::endl;
   std::cout << "Example:" << std::endl;
   std::cout << "   ./gridpp_train testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=parameters.txt -c regression -v T" << std::endl;
   std::cout << std::endl;
   std::cout << "Obs/fcst files" << std::endl;
   std::cout << "   Obs/fcst file types are autodetected, but can be specified using:" << std::endl;
   std::cout << File::getDescriptions();
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
   double startTime = Util::clock();
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
   int nTime = forecast->getNumTime();
   vec2 lats = forecast->getLats();
   vec2 lons = forecast->getLons();
   vec2 elevs = forecast->getElevs();

   Variable::Type variable = Variable::T;
   std::vector<double> offsets = forecast->getTimes();
   int D = setup.forecasts.size();

   std::cout << "Number of obs files: " << setup.observations.size() << std::endl;
   std::cout << "Number of fcst files: " << setup.forecasts.size() << std::endl;
   std::cout << "Number of lats: " << nLat << std::endl;
   std::cout << "Number of lons: " << nLon << std::endl;
   std::cout << "Number of times: " << nTime << std::endl;
   std::cout << "Forecast reference times: ";
   for(int f = 0; f < setup.forecasts.size(); f++) {
      double timeF = setup.forecasts[f]->getReferenceTime();
      std::cout << " " << timeF;
   }
   std::cout << std::endl;
   std::cout << "Observation reference times: ";
   for(int f = 0; f < setup.observations.size(); f++) {
      double timeO = setup.observations[f]->getReferenceTime();
      std::cout << " " << timeO;
   }
   std::cout << std::endl;

   // figure out which files match
   std::vector<int> indexF;
   std::vector<int> indexO;
   for(int f = 0; f < setup.forecasts.size(); f++) {
      int io = Util::MV;
      double timeF = setup.forecasts[f]->getReferenceTime();
      for(int o = 0; o < setup.observations.size(); o++) {
         double timeO = setup.observations[o]->getReferenceTime();
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
   std::vector<double> ftimes = forecast->getTimes();
   std::vector<double> otimes = observation->getTimes();
   double freftime = forecast->getReferenceTime();
   double oreftime = observation->getReferenceTime();

   std::cout << "Forecast offsets:";
   for(int t = 0; t < ftimes.size(); t++) {
      std::cout << " " << ftimes[t] - freftime;
   }
   std::cout << std::endl;
   std::cout << "Observation offsets:";
   for(int t = 0; t < otimes.size(); t++) {
      std::cout << " " << otimes[t] - oreftime;
   }
   std::cout << std::endl;

   // If multiple available obstimes, then use the one with the lowest ftime
   for(int t = 0; t < offsets.size(); t++) {
      std::cout << "Forecast time index: " << t << std::endl;
      double ftime = ftimes[t] - freftime;
      int fTimeIndex = t;
      int oTimeIndex = Util::MV;
      int fmaxDays = ftimes[ftimes.size()-1] / 86400;
      for(int d = 0; d < fmaxDays; d++) {
         // Find appropriate obs time
         for(int tt = 0; tt < otimes.size(); tt++) {
            double otime = otimes[tt] - oreftime + d * 86400;
            if(otime == ftime) {
               oTimeIndex = tt;
               break;
            }
         }
      }
      assert(Util::isValid(oTimeIndex));

      // Loop over all files
      std::vector<FieldPtr> ffields;
      std::vector<FieldPtr> ofields;
      for(int d = 0; d < setup.forecasts.size(); d++) {
         double ftime = setup.forecasts[d]->getTimes()[t];
         int fDateIndex = d;
         int oDateIndex = Util::MV;

         // Find corresponding observation file
         for(int dd = 0; dd < setup.observations.size(); dd++) {
            double otime = setup.observations[dd]->getTimes()[oTimeIndex];
            if(ftime == otime) {
               oDateIndex = dd;
               break;
            }
         }

         if(Util::isValid(fDateIndex) && Util::isValid(oDateIndex)) {
            std::cout << "   Found matching: (" << fDateIndex << "," << fTimeIndex << ") " 
                      << " (" << oDateIndex << "," << oTimeIndex << ")" << std::endl;
            assert(setup.forecasts.size() > fDateIndex);
            assert(setup.observations.size() > oDateIndex);
            ffields.push_back(setup.forecasts[fDateIndex]->getField(variable, fTimeIndex));
            ofields.push_back(setup.observations[oDateIndex]->getField(variable, oTimeIndex));
         }
      }
      if(ffields.size() == 0) {
         std::stringstream ss;
         ss << "Cannot find files with available data for timestep " << t << std::endl;
         Util::error(ss.str());
      }
      assert(ffields.size() > 0);

      ////////////////////////
      // Compute parameters
      // Do this in parallel, inserting the values into a preallocated matrix
      std::vector<std::vector<std::vector<float> > > parameters;
      parameters.resize(nLat);
      for(int i = 0; i < nLat; i++) {
         parameters[i].resize(nLon);
      }

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            // Arrange data
            std::vector<ObsEns> data;
            for(int d = 0; d < ffields.size(); d++){
               // Downscaling (currently nearest neighbour)
               float obs = (*ofields[d])(I[i][j],J[i][j],0);
               const Ens& ens = (*ffields[d])(i,j);
               ObsEns obsens(obs, ens);
               data.push_back(obsens);
            }

            Parameters par = setup.method->train(data);
            parameters[i][j] = par.getValues();
         }
      }

      // Then save values serially
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int offset = t;
            float lat = lats[i][j];
            float lon = lons[i][j];
            float elev = elevs[i][j];

            setup.output->setParameters(Parameters(parameters[i][j]), offset, Location(lat, lon, elev));
         }
      }
   }

   std::cout << "Writing to parameter file" << std::endl;
   setup.output->write();
   double endTime = Util::clock();
   std::cout << "Total time: " << endTime - startTime << std::endl;
}
