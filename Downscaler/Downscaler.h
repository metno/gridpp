#ifndef DOWNSCALER_H
#define DOWNSCALER_H
#include "../DataFile.h"
#include "../ParameterFile.h"

class Downscaler {
   public:
      Downscaler(const ParameterFile& iParameterFile);
      void downscale(const DataFile& iInput, DataFile& iOutput) const;
   protected:
      void downscaleCore(const DataFile& iInput, DataFile& iOutput) const;
};
#endif
