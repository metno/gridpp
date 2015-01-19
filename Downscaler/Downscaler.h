#ifndef DOWNSCALER_H
#define DOWNSCALER_H
class File;

class Downscaler {
   public:
      Downscaler();
      void downscale(const File& iInput, File& iOutput) const;
   protected:
      virtual void downscaleCore(const File& iInput, File& iOutput) const = 0;
};
#include "NearestNeighbour.h"
#include "Gradient.h"
#include "Smart.h"
#endif
