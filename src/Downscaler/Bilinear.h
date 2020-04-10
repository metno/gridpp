#ifndef DOWNSCALER_BILINEAR_H
#define DOWNSCALER_BILINEAR_H
#include <map>
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
class File;
class DownscalerBilinear : public Downscaler {
   public:
      DownscalerBilinear(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "bilinear";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
};
#endif
