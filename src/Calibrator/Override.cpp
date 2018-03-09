#include "Override.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorOverride::CalibratorOverride(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mRadius(0) {
   iOptions.getValue("radius", mRadius);
   iOptions.check();
}
bool CalibratorOverride::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   int nX = iFile.getNumX();
   int nY = iFile.getNumY();
   vec2 elevs = iFile.getElevs();

   if(iParameterFile->getNumParameters() != 1) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must have two dataacolumns");
   }
   if(!iParameterFile->isLocationDependent()) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must be location dependent");
   }

   std::vector<Location> pointLocations = iParameterFile->getLocations();
   KDTree searchTree(iFile.getLats(), iFile.getLons());
   std::vector<int> Is(pointLocations.size(), 0);
   std::vector<int> Js(pointLocations.size(), 0);
   for(int k = 0; k < pointLocations.size(); k++) {
      searchTree.getNearestNeighbour(pointLocations[k].lat(), pointLocations[k].lon(), Is[k], Js[k]);
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      const FieldPtr field = iFile.getField(mVariable, t);
      for(int k = 0; k < pointLocations.size(); k++) {
         float currElev = pointLocations[k].elev();
         Parameters parameters = iParameterFile->getParameters(t, pointLocations[k]);
         float value = parameters[0];
         int Inn = Is[k];
         int Jnn = Js[k];
         for(int e = 0; e < nEns; e++) {
            if(Util::isValid(currElev) && mRadius > 0 && Inn >= mRadius && Inn < nX - mRadius && Jnn >= mRadius && Jnn < nY - mRadius) {
               // Search in a neighbourhood
               int Ibest = Inn;
               int Jbest = Jnn;
               float minElevDiff = Util::MV;
               for(int I = Inn - mRadius; I <= Inn + mRadius; I++) {
                  for(int J = Jnn - mRadius; J <= Jnn + mRadius; J++) {
                     if(Util::isValid(elevs[I][J])) {
                        float elevDiff = abs(elevs[I][J] - currElev);
                        if(!Util::isValid(minElevDiff) || elevDiff < minElevDiff) {
                           Ibest = I;
                           Jbest = J;
                           minElevDiff = elevDiff;
                        }
                     }
                  }
               }
               (*field)(Ibest, Jbest, e) = value;
            }
            else {
               (*field)(Inn, Jnn, e) = value;
            }
         }
      }
   }
   return true;
}

std::string CalibratorOverride::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c override", "Overrides certain points in the file based on parameter values at point locations.") << std::endl;
   ss << Util::formatDescription("   radius=0", "Look in a square box with this radius (in number of gridpoints) around the point and override the gridpoint with the closest elevation to the elevation in the parameter file.") << std::endl;
   return ss.str();
}
