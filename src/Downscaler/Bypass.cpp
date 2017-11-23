#include "Bypass.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerBypass::DownscalerBypass(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
}

void DownscalerBypass::downscaleCore(const File& iInput, File& iOutput) const {
   if(!iOutput.hasVariable(mOutputVariable)) {
      // int nTime = iInput.getNumTime();
      // for(int t = 0; t < nTime; t++) {
      iOutput.initNewVariable(mOutputVariable);
      // }
   }
}

std::string DownscalerBypass::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-d bypass", "No downscaling is done, but useful if the variable is derived by a calibrator.") << std::endl;
   return ss.str();
}
