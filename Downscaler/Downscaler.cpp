#include "Downscaler.h"

Downscaler::Downscaler(Variable::Type iVariable) :
      mVariable(iVariable) {
}

void Downscaler::downscale(const File& iInput, File& iOutput) const {
   downscaleCore(iInput, iOutput);
}

Downscaler* Downscaler::getScheme(std::string iType, Variable::Type iVariable, std::string iOptions) {
   if(iType == "nearestneighbour") {
      return new DownscalerNearestNeighbour(iVariable);
   }
   else if(iType == "gradient") {
      return new DownscalerGradient(iVariable);
   }
   else if(iType == "smart") {
      return new DownscalerSmart(iVariable);
   }
   else {
      Util::error("Could not instantiate downscaler of type '" + iType + "'");
      return NULL;
   }
}
