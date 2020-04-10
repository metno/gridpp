#include "Bilinear.h"
#include "../File/File.h"
#include "../Util.h"
#include "gridpp.h"
#include <math.h>

DownscalerBilinear::DownscalerBilinear(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
   iOptions.check();
}

void DownscalerBilinear::downscaleCore(const File& iInput, File& iOutput) const {
   int nEns = iInput.getNumEns();
   int nTime = iInput.getNumTime();

   gridpp::Grid igrid(iInput.getLats(), iInput.getLons());
   gridpp::Grid ogrid(iOutput.getLats(), iOutput.getLons());
   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);
      for(int e = 0; e < nEns; e++) {
         ofield.set(gridpp::bilinear(igrid, ogrid, ifield(e)), e);
      }
   }
}

std::string DownscalerBilinear::description(bool full) {
   std::stringstream ss;
   if(full)
      ss << Util::formatDescription("-d bilinear", "Bilinear interpolation using the 4 points surrounding the lookup point. If the lookup point is outside the input domain, then the nearest neighbour is used.") << std::endl;
   else
      ss << Util::formatDescription("-d bilinear", "Bilinear interpolation") << std::endl;
   return ss.str();
}
