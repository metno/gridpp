#include "Mask.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
CalibratorMask::CalibratorMask(Variable::Type iVariable, const Options& iOptions) :
      Calibrator(iOptions),
      mVariable(iVariable),
      mUseNearestOnly(false) {
}
bool CalibratorMask::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(!iParameterFile->isLocationDependent()) {
      Util::error("Cannot use a location independent parameter file in mask calibrator");
   }
   if(!iParameterFile->isFixedSize()) {
      Util::error("Cannot use a parameter file without a constant number of parameters");
   }
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();
   vec2 keep;
   keep.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      // Default to not keep, unless we find a parameter point within the radius
      keep[i].resize(nLon, 0);
   }

   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            if(t == 0 || iParameterFile->isTimeDependent()) {
               Location loc(Util::MV, Util::MV, Util::MV);
               Location currLocation(lats[i][j], lons[i][j], elevs[i][j]);
               if(mUseNearestOnly) {
                  iParameterFile->getNearestLocation(t, currLocation, loc);
                  Parameters parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
                  float dist = currLocation.getDistance(loc);
                  if(Util::isValid(dist) && dist < parameters[0]) {
                     keep[i][j] = 1;
                  }
               }
               else {
                  std::vector<Location> locations = iParameterFile->getLocations();
                  for(int k = 0; k < locations.size(); k++) {
                     Parameters parameters = iParameterFile->getParameters(t, locations[k]);
                     float dist = currLocation.getDistance(locations[k]);
                     if(Util::isValid(dist) && dist < parameters[0]) {
                        keep[i][j] = 1;
                     }
                  }
               }
            }
            if(keep[i][j] == 0) {
               for(int e = 0; e < nEns; e++) {
                  (*field)(i,j,e) = Util::MV;
               }
            }
         }
      }
   }
   iFile.setElevs(elevs);
   return true;
}

std::string CalibratorMask::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c mask", "Only keep gridpoints that are within a certain radius of points in the parameter file. All other gridpoints are set to missing. The parameter file must be spatial and conain one parameter representing the radius in meters. Each gridpoint is checked against each parameter point to see if the gridpoint is within its radius.") << std::endl;
   return ss.str();
}
