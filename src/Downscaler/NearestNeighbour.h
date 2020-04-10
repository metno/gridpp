#ifndef DOWNSCALER_NEAREST_NEIGHBOUR_H
#define DOWNSCALER_NEAREST_NEIGHBOUR_H
#include <map>
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
class File;
class DownscalerNearestNeighbour : public Downscaler {
   public:
      DownscalerNearestNeighbour(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "nearestNeighbour";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
};
#endif
