#include "Smart.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerSmart::DownscalerSmart(Variable::Type iVariable) :
      Downscaler(iVariable),
      mSearchRadius(3),
      mNumSmart(5) {
}

void DownscalerSmart::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumLat();
   int nLon = iOutput.getNumLon();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   vec2 ilats  = iInput.getLats();
   vec2 ilons  = iInput.getLons();
   vec2 ielevs = iInput.getElevs();
   vec2 olats  = iOutput.getLats();
   vec2 olons  = iOutput.getLons();
   vec2 oelevs = iOutput.getElevs();

   // Get nearest neighbour
   vec3Int nearestI, nearestJ;
   getSmartNeighbours(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mVariable, t);
      Field& ofield = *iOutput.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               float total = 0;
               int   count = 0;
               for(int n = 0; n < mNumSmart; n++) {
                  int ii = nearestI[i][j][n];
                  int jj = nearestJ[i][j][n];

                  if(Util::isValid(ii) && Util::isValid(jj)) {
                     float value = ifield[ii][jj][e];
                     if(Util::isValid(value)) {
                        total += value;
                        count++;
                     }
                  }
               }
               if(count > 0)
                  ofield[i][j][e] = total/count;
               else
                  ofield[i][j][e] = Util::MV;
               // ofield[i][j][e] = nearestI[i][j][0];
            }
         }
      }
   }
}
void DownscalerSmart::setNumSmart(int iNumSmart) {
   mNumSmart = iNumSmart;
}
void DownscalerSmart::setSearchRadius(int iNumPoints) {
   mSearchRadius = iNumPoints;
}
void DownscalerSmart::getSmartNeighbours(const File& iFrom, const File& iTo, vec3Int& iI, vec3Int& iJ) const {
   getSmartNeighbours(iFrom, iTo, mSearchRadius, mNumSmart, iI, iJ);
}

void DownscalerSmart::getSmartNeighbours(const File& iFrom, const File& iTo, int iSearchRadius, int iNumSmart, vec3Int& iI, vec3Int& iJ) {
   vec2 ilats  = iFrom.getLats();
   vec2 ilons  = iFrom.getLons();
   vec2 ielevs = iFrom.getElevs();
   vec2 olats  = iTo.getLats();
   vec2 olons  = iTo.getLons();
   vec2 oelevs = iTo.getElevs();
   int nLon    = iTo.getNumLon();
   int nLat    = iTo.getNumLat();
   int numSearch = getNumSearchPoints(iSearchRadius);

   vec2Int Icenter, Jcenter;
   DownscalerNearestNeighbour::getNearestNeighbourFast(iFrom, iTo, Icenter, Jcenter);

   iI.resize(nLat);
   iJ.resize(nLat);
   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      iI[i].resize(nLon);
      iJ[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         int Ic = Icenter[i][j];
         int Jc = Jcenter[i][j];
         float oelev = ielevs[i][j];

         // Compute elevation differences on the stencil surrounding the current point
         std::vector<std::pair<int, float> > elevDiff; // elevation difference of points in stencil
                                                       // but placed in 1D array, to simplify sort
         elevDiff.reserve(numSearch);
         // Keep track of which I/J corresponds to indices in elevDiff
         std::vector<int> Ilookup;
         std::vector<int> Jlookup;
         Ilookup.reserve(numSearch);
         Jlookup.reserve(numSearch);

         int index = 0;
         for(int ii = std::max(0, Ic-iSearchRadius); ii <= std::min(iFrom.getNumLat()-1, Ic+iSearchRadius); ii++) {
            for(int jj = std::max(0, Jc-iSearchRadius); jj <= std::min(iFrom.getNumLon()-1, Jc+iSearchRadius); jj++) {
               float ielev = ielevs[ii][jj];
               float diff = 1e10;
               if(Util::isValid(ielev) && Util::isValid(oelev))
                  diff = abs(ielev - oelev);
               elevDiff.push_back(std::pair<int, float>(index, diff));
               Ilookup.push_back(ii);
               Jlookup.push_back(jj);
               index++;
            }
         }

         std::sort(elevDiff.begin(), elevDiff.end(), Util::sort_pair_second<int,float>());

         std::stringstream ss;
         ss << "Smart neighbours for " << i << " " << j << " " << oelev;
         Util::status(ss.str());

         iI[i][j].resize(iNumSmart, Util::MV);
         iJ[i][j].resize(iNumSmart, Util::MV);
         // Use nearest neighbour if all fails
         if(elevDiff.size() == 0) {
            iI[i][j][0] = Ic;
            iJ[i][j][0] = Jc;
         }
         else {
            for(int n = 0; n < iNumSmart; n++) {
               if(elevDiff.size() > n) {
                  int index = elevDiff[n].first;
                  iI[i][j][n] = Ilookup[index];
                  iJ[i][j][n] = Jlookup[index];
                  std::stringstream ss;
                  ss << "   " << iI[i][j][n] << " " << iJ[i][j][n] << " " << elevDiff[n].second;
                  Util::status(ss.str());
               }
            }
         }
      }
   }
}

int DownscalerSmart::getNumSearchPoints() const {
   return getNumSearchPoints(mSearchRadius);
}

int DownscalerSmart::getNumSearchPoints(int iSearchRadius) {
   return (iSearchRadius*2+1)*(iSearchRadius*2+1);
}

int  DownscalerSmart::getSearchRadius() const {
   return mSearchRadius;
}

std::string DownscalerSmart::description() {
   std::stringstream ss;
   ss << "   -d smart" << std::endl;
   ss << "      searchRadius=3            Search for smart neighbours within this radius (gridpoints)" << std::endl;
   ss << "      numSmart=5                Average this many smart neighbours" << std::endl;
   return ss.str();
}
