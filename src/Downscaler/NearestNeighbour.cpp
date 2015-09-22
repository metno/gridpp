#include "NearestNeighbour.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

// std::map<const File*, std::map<const File*, std::pair<vec2Int, vec2Int> > > DownscalerNearestNeighbour::mNeighbourCache;

DownscalerNearestNeighbour::DownscalerNearestNeighbour(Variable::Type iVariable, const Options& iOptions) :
      Downscaler(iVariable, iOptions) {
}

void DownscalerNearestNeighbour::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumLat();
   int nLon = iOutput.getNumLon();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   getNearestNeighbour(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mVariable, t);
      Field& ofield = *iOutput.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int I = nearestI[i][j];
            int J = nearestJ[i][j];
            for(int e = 0; e < nEns; e++) {
               if(Util::isValid(I) && Util::isValid(J))
                  ofield(i,j,e) = ifield(I,J,e);
               else
                  ofield(i,j,e) = Util::MV;
            }
         }
      }
   }
}

std::string DownscalerNearestNeighbour::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-d nearestNeighbour", "Uses the nearest gridpoint in curved distance") << std::endl;
   return ss.str();
}
