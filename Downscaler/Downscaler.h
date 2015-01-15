#ifndef DOWNSCALER_H
#define DOWNSCALER_H
#include "../File.h"
#include "../ParameterFile.h"

class Downscaler {
   public:
      Downscaler(const ParameterFile& iParameterFile);
      void downscale(const File& iInput, File& iOutput) const;
   protected:
      void downscaleCore(const File& iInput, File& iOutput) const;
};
#endif
