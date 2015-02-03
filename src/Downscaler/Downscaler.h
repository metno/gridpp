#ifndef DOWNSCALER_H
#define DOWNSCALER_H
#include <string>
#include <map>
#include <boost/uuid/uuid.hpp>
#include "../Options.h"
#include "../Variable.h"
class File;
typedef std::vector<std::vector<int> > vec2Int;

class Downscaler {
   public:
      Downscaler(Variable::Type iVariable);
      virtual ~Downscaler() {};
      bool downscale(const File& iInput, File& iOutput) const;
      static Downscaler* getScheme(std::string iName, Variable::Type iVariable, const Options& iOptions);
      virtual std::string name() const = 0;

      // Slow method: Check every combination
      // Return Util::MV when it cannot find a neighbour
      static void getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);
      // Faster method: Assume lats/lons are sorted
      static void getNearestNeighbourFast(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);
   protected:
      virtual void downscaleCore(const File& iInput, File& iOutput) const = 0;
      Variable::Type mVariable;
   private:
      // Cache calls to nearest neighbour
      //! Is the nearest neighbours in @param iFrom for each point in @param iTo already computed?
      static bool isCached(const File& iFrom, const File& iTo);
      static void addToCache(const File& iFrom, const File& iTo, vec2Int iI, vec2Int iJ);
      static bool getFromCache(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);
      static std::map<boost::uuids::uuid, std::map<boost::uuids::uuid, std::pair<vec2Int, vec2Int> > > mNeighbourCache;
};
#include "NearestNeighbour.h"
#include "Gradient.h"
#include "Smart.h"
#include "Bypass.h"
#include "Pressure.h"
#endif
