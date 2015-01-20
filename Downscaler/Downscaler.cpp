#include "Downscaler.h"

Downscaler::Downscaler(Variable::Type iVariable) :
      mVariable(iVariable) {
}

void Downscaler::downscale(const File& iInput, File& iOutput) const {
   downscaleCore(iInput, iOutput);
}

Downscaler* Downscaler::getScheme(std::string iType, Variable::Type iVariable, Options& iOptions) {
   if(iType == "nearestneighbour") {
      return new DownscalerNearestNeighbour(iVariable);
   }
   else if(iType == "gradient") {
      DownscalerGradient* d = new DownscalerGradient(iVariable);
      float constantGradient = Util::MV;
      if(iOptions.getValue("constantGradient", constantGradient)) {
         d->setConstantGradient(constantGradient);
      }
      float searchRadius = Util::MV;
      if(iOptions.getValue("searchRadius", searchRadius)) {
         d->setSearchRadius(searchRadius);
      }
      return d;
   }
   else if(iType == "smart") {
      return new DownscalerSmart(iVariable);
   }
   else {
      Util::error("Could not instantiate downscaler of type '" + iType + "'");
      return NULL;
   }
}
