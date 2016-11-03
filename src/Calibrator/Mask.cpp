#include "Mask.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
CalibratorMask::CalibratorMask(Variable::Type iVariable, const Options& iOptions) :
      Calibrator(iOptions), mVariable(iVariable) {
}
bool CalibratorMask::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(!iParameterFile->isLocationDependent()) {
      Util::error("Cannot use a location independent parameter file in mask calibrator");
   }
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();
   vec2 inMask;
   inMask.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      inMask[i].resize(nLon, 0);
   }

   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            if(t == 0 || iParameterFile->isTimeDependent()) {
               Location loc(Util::MV, Util::MV, Util::MV);
               Location currLocation(lats[i][j], lons[i][j], elevs[i][j]);
               iParameterFile->getNearestLocation(t, currLocation, loc);
               Parameters parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
               float dist = currLocation.getDistance(loc);
               if(Util::isValid(dist) && dist > parameters[0]) {
                  inMask[i][j] = 1;
               }
            }
            if(inMask[i][j] == 1) {
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
   ss << Util::formatDescription("-c altitude", "Changes the altitudes in the file to the altitudes in the parameter file. If the file does not have an altitude field, then a new one is not created.") << std::endl;
   return ss.str();
}
