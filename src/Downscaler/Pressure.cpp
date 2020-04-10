#include "Pressure.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>
const float DownscalerPressure::mConstant = -1.21e-4;

DownscalerPressure::DownscalerPressure(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
   iOptions.check();
}

void DownscalerPressure::downscaleCore(const File& iInput, File& iOutput) const {
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   gridpp::Grid igrid(iInput.getLats(), iInput.getLons());
   gridpp::Grid ogrid(iOutput.getLats(), iOutput.getLons());
   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);

      for(int e = 0; e < nEns; e++) {
         ofield.set(gridpp::pressure(igrid, ogrid, ifield(e)), e);
      }
   }
}
std::string DownscalerPressure::description(bool full) {
   std::stringstream ss;
   if(full)
      ss << Util::formatDescription("-d pressure", "Adjusts the pressure of the nearest neighbour based on the elevation difference and a standard atmosphere.") << std::endl;
   else
      ss << Util::formatDescription("-d pressure", "Adjusts the pressure based on the elevation differences") << std::endl;
   return ss.str();
}

float DownscalerPressure::calcPressure(float iElev0, float iPressure0, float iElev1) {
   if(Util::isValid(iElev0) && Util::isValid(iPressure0) && Util::isValid(iElev1)) {
      float dElev = iElev1 - iElev0;
      return iPressure0 * exp(mConstant * (dElev));
   }
   else {
      return Util::MV;
   }
}
