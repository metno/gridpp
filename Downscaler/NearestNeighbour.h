#ifndef DOWNSCALER_NEAREST_NEIGHBOUR_H
#define DOWNSCALER_NEAREST_NEIGHBOUR_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
class DownscalerNearestNeighbour : public Downscaler {
   public:
      DownscalerNearestNeighbour(Variable::Type iVariable);
      // Slow method: Check every combination
      static void getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);
      // Faster method: Assume lats/lons are sorted
      static void getNearestNeighbourFast(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);
      static std::string description();
      std::string name() const {return "nearestNeighbour";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
};
#endif
