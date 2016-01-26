#include "Altitude.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
CalibratorAltitude::CalibratorAltitude(const Options& iOptions) :
      Calibrator(iOptions) {
}
bool CalibratorAltitude::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(!iParameterFile->isLocationDependent()) {
      Util::error("Cannot use a location independent parameter file to update the altitudes");
   }
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         Location loc(Util::MV, Util::MV, Util::MV);
         iParameterFile->getNearestLocation(0, Location(lats[i][j], lons[i][j], elevs[i][j]), loc);
         elevs[i][j] = loc.elev();
      }
   }
   iFile.setElevs(elevs);
   return true;
}

std::string CalibratorAltitude::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c altitude", "Changes the altitudes in the file to the altitudes in the parameter file.") << std::endl;
   return ss.str();
}
