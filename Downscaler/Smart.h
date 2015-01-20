#ifndef DOWNSCALER_SMART_H
#define DOWNSCALER_SMART_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
typedef std::vector<std::vector<std::vector<int> > > vec3Int;
class DownscalerSmart : public Downscaler {
   public:
      //! Downscale the specified variable
      DownscalerSmart(Variable::Type iVariable);
      //! Use this many smart neighbours
      void setNumSmart(int iNumSmart);
      //! Search for smart neighbours within a neighbourhood of points within +- iNumPoints
      //! in each direction.
      void setSearchRadius(int iNumPoints);

      static void getSmartNeighbours(const File& iFrom, const File& iTo, int iSearchRadius, int iNumSmart, vec3Int& iI, vec3Int& iJ);
             void getSmartNeighbours(const File& iFrom, const File& iTo, vec3Int& iI, vec3Int& iJ) const;
      static int getNumSearchPoints(int iSearchRadius) ;
             int getNumSearchPoints() const;
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      int mSearchRadius;
      int mNumSmart;
};
#endif
