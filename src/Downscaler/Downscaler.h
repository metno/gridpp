#ifndef DOWNSCALER_H
#define DOWNSCALER_H
#include <string>
#include <map>
#include <boost/uuid/uuid.hpp>
#include "../Options.h"
#include "../Variable.h"
class File;
typedef std::vector<std::vector<int> > vec2Int;

//! Converts fields from one grid to another
class Downscaler {
   public:
      Downscaler(Variable::Type iVariable);
      virtual ~Downscaler() {};
      bool downscale(const File& iInput, File& iOutput) const;
      static Downscaler* getScheme(std::string iName, Variable::Type iVariable, const Options& iOptions);
      virtual std::string name() const = 0;

      //! Create a nearest-neighbour map. For each grid point in iTo, find the index into the grid
      //! in iFrom of the nearest neighbour. Uses a brute force method by checking every
      //! neighbour.
      //! @param iI I-indices of nearest point. Set to Util::MV if no nearest neighbour.
      //! @param iJ J-indices of nearest point. Set to Util::MV if no nearest neighbour.
      static void getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);

      //! A faster version of getNearestNeighbour, which assuming that latitudes and longitudes in
      //! iFrom are sorted. Defaults to the slow method if lats/lons are not sorted.
      static void getNearestNeighbourFast(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ);

      //! Find the index into the lat/lon arrays in iFrom that is nearest to the point
      //! defined by iLon,iLat
      //! @param iI I-index of nearest point. Set to Util::MV if no nearest neighbour.
      //! @param iJ J-index of nearest point. Set to Util::MV if no nearest neighbour.
      static void getNearestNeighbour(const File& iFrom, float iLon, float iLat, int& iI, int &iJ);
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
