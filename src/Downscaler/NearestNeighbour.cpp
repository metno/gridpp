#include "NearestNeighbour.h"
#include "../File/File.h"
#include "../Util.h"
#include "gridpp.h"
#include <math.h>

DownscalerNearestNeighbour::DownscalerNearestNeighbour(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
   iOptions.check();
}

void DownscalerNearestNeighbour::downscaleCore(const File& iInput, File& iOutput) const {
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   gridpp::Grid igrid(iInput.getLats(), iInput.getLons());
   gridpp::Grid ogrid(iOutput.getLats(), iOutput.getLons());
   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);

      for(int e = 0; e < nEns; e++) {
         ofield.set(gridpp::nearest(igrid, ogrid, ifield(e)), e);
      }
   }
}
std::string DownscalerNearestNeighbour::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-d nearestNeighbour", "Uses the nearest gridpoint in curved distance") << std::endl;
   return ss.str();
}
