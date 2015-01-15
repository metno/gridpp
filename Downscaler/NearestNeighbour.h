#ifndef DOWNSCALER_NEAREST_NEIGHBOUR_H
#define DOWNSCALER_NEAREST_NEIGHBOUR_H
#include "Downscaler.h"

class DownscalerNearestNeighbour : public Downscaler {
   public:
      DownscalerNearestNeighbour(const ParameterFile& iParameterFile);
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
};
#endif
