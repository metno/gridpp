#include "Override.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorOverride::CalibratorOverride(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mRadius(0),
      mMaxElevDiff(Util::MV) {
   iOptions.getValue("radius", mRadius);
   iOptions.getValue("maxElevDiff", mMaxElevDiff);
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
            for(int I = std::max(0, Inn - mRadius); I <= std::min(nX-1, Inn + mRadius); I++) {
               for(int J = std::max(0, Jnn - mRadius); J <= std::min(nY-1, Jnn + mRadius); J++) {
                  if(!Util::isValid(mMaxElevDiff)) {
                     // Ignore elevation information
                     (*field)(I, J, e) = value;
                  }
                  else if(Util::isValid(elevs[I][J]) && Util::isValid(currElev)) {
                     // Check if within elevation bound
                     float elevDiff = abs(elevs[I][J] - currElev);
                     if(elevDiff < mMaxElevDiff) {
                        (*field)(I, J, e) = value;
                     }
                  }
                  else {
                     // Do not update
                  }
               }
            }
         }
      }
   }
   return true;
}

std::string CalibratorOverride::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c override", "Overrides certain points in the file based on parameter values at point locations.") << std::endl;
   ss << Util::formatDescription("   radius=0", "Write values into all gridpoints in a square box with this radius (in number of gridpoints).") << std::endl;
   ss << Util::formatDescription("   maxElevDiff=undef", "Only write in gridpoints that are within this elevation difference to the point. If undefined, then do not do an elevation check.") << std::endl;
   return ss.str();
}
