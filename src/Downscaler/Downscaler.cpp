#include "Downscaler.h"
#include <cmath>
#include "../File/File.h"

std::map<boost::uuids::uuid, std::map<boost::uuids::uuid, std::pair<vec2Int, vec2Int> > > Downscaler::mNeighbourCache;

Downscaler::Downscaler(Variable::Type iVariable) :
      mVariable(iVariable) {
}

bool Downscaler::downscale(const File& iInput, File& iOutput) const {
   if(iInput.getNumTime() != iOutput.getNumTime())
      return false;

   downscaleCore(iInput, iOutput);
   return true;
}

Downscaler* Downscaler::getScheme(std::string iName, Variable::Type iVariable, const Options& iOptions) {
   if(iName == "nearestNeighbour") {
      return new DownscalerNearestNeighbour(iVariable);
   }
   else if(iName == "gradient") {
      DownscalerGradient* d = new DownscalerGradient(iVariable);
      float constantGradient = Util::MV;
      if(iOptions.getValue("constantGradient", constantGradient)) {
         d->setConstantGradient(constantGradient);
      }
      float searchRadius = Util::MV;
      if(iOptions.getValue("searchRadius", searchRadius)) {
         d->setSearchRadius(searchRadius);
      }
      return d;
   }
   else if(iName == "smart") {
      DownscalerSmart* d = new DownscalerSmart(iVariable);
      float searchRadius = Util::MV;
      if(iOptions.getValue("searchRadius", searchRadius)) {
         d->setSearchRadius(searchRadius);
      }
      float numSmart = Util::MV;
      if(iOptions.getValue("numSmart", numSmart)) {
         d->setNumSmart(numSmart);
      }
      return d;
   }
   else if(iName == "bypass") {
      DownscalerBypass* d = new DownscalerBypass(iVariable);
      return d;
   }
   else {
      Util::error("Could not instantiate downscaler of type '" + iName + "'");
      return NULL;
   }
}

void Downscaler::getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(isCached(iFrom, iTo)) {
       getFromCache(iFrom, iTo, iI, iJ);
       return;
   }

   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumLon();
   int nLat = iTo.getNumLat();

   iI.resize(nLat);
   iJ.resize(nLat);

   // Check if the grid is the same
   if(iFrom.getNumLat() == iTo.getNumLat() && iFrom.getNumLon() == iTo.getNumLon()) {
      if(ilats == olats && ilons == olons) {
         for(int i = 0; i < nLat; i++) {
            iI[i].resize(nLon, 0);
            iJ[i].resize(nLon, 0);
            for(int j = 0; j < nLon; j++) {
               iI[i][j] = i;
               iJ[i][j] = j;
            }
         }
      }
      Util::status("Grids are identical, short cut in finding nearest neighbours");
      addToCache(iFrom, iTo, iI, iJ);
      return;
   }

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      iI[i].resize(nLon, Util::MV);
      iJ[i].resize(nLon, Util::MV);
      for(int j = 0; j < nLon; j++) {
         if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
            float minDist = Util::MV;
            int I = Util::MV;
            int J = Util::MV;
            for(int ii = 0; ii < iFrom.getNumLat(); ii++) {
               for(int jj = 0; jj < iFrom.getNumLon(); jj++) {
                  if(Util::isValid(ilats[ii][jj]) && Util::isValid(ilons[ii][jj])) {
                     float currDist = Util::getDistance(olats[i][j], olons[i][j], ilats[ii][jj], ilons[ii][jj]);
                     if(!Util::isValid(minDist) || currDist < minDist) {
                        I = ii;
                        J = jj;
                        minDist = currDist;
                     }
                  }
               }
            }
            iI[i][j] = I;
            iJ[i][j] = J;
         }
      }
   }
   addToCache(iFrom, iTo, iI, iJ);
}

void Downscaler::getNearestNeighbourFast(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(iTo.getNumLat() == 0 || iTo.getNumLon() == 0) {
      return;
   }

   if(isCached(iFrom, iTo)) {
      getFromCache(iFrom, iTo, iI, iJ);
      return;
   }

   if(iFrom.getNumLat() == 0 || iFrom.getNumLon() == 0) {
      iI.resize(iTo.getNumLat());
      iJ.resize(iTo.getNumLat());
      for(int i = 0; i < iTo.getNumLat(); i++) {
         iI[i].resize(iTo.getNumLon(), Util::MV);
         iJ[i].resize(iTo.getNumLon(), Util::MV);
      }
      return;
   }

   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumLon();
   int nLat = iTo.getNumLat();

   // Check if grid has missing points
   bool hasMissing = false;
   for(int ii = 0; ii < iFrom.getNumLat(); ii++) {
      for(int jj = 0; jj < iFrom.getNumLon(); jj++) {
         if(!Util::isValid(ilats[ii][jj]) || !Util::isValid(ilons[ii][jj])) {
            hasMissing = true;
         }
      }
   }
   if(hasMissing) {
      std::stringstream ss;
      ss << "Lats/lons has missing gridpoints, using a slower method to find nearest neighbours";
      Util::warning(ss.str());
      return getNearestNeighbour(iFrom, iTo, iI, iJ);
   }
   // Check if grid is sorted
   bool isSorted = true;
   for(int ii = 1; ii < iFrom.getNumLat(); ii++) {
      for(int jj = 1; jj < iFrom.getNumLon(); jj++) {
         if(!Util::isValid(ilats[ii][jj]) || !Util::isValid(ilons[ii][jj]) || ilats[ii][jj] < ilats[ii-1][jj] || ilons[ii][jj] < ilons[ii][jj-1]) {
            isSorted = false;
         }
      }
   }
   if(!isSorted) {
      std::stringstream ss;
      ss << "Lats/lons are not sorted, using a slower method to find nearest neighbours";
      Util::warning(ss.str());
      return getNearestNeighbour(iFrom, iTo, iI, iJ);
   }

   iI.resize(nLat);
   iJ.resize(nLat);

   // Check if the grid is the same
   if(iFrom.getNumLat() == iTo.getNumLat() && iFrom.getNumLon() == iTo.getNumLon()) {
      if(ilats == olats && ilons == olons) {
         for(int i = 0; i < nLat; i++) {
            iI[i].resize(nLon, 0);
            iJ[i].resize(nLon, 0);
            for(int j = 0; j < nLon; j++) {
               iI[i][j] = i;
               iJ[i][j] = j;
            }
         }
      }
      Util::status("Grids are identical, short cut in finding nearest neighbours");
      return;
   }

   float tol = 0.2;

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      int I = 0;
      int J = 0;
      iI[i].resize(nLon, 0);
      iJ[i].resize(nLon, 0);
      for(int j = 0; j < nLon; j++) {
         // if(j % 10 == 0)
         //    std::cout << i << " " << j << std::endl;
         int counter = 0;
         float currLat = olats[i][j];
         float currLon = olons[i][j];
         while(true) {
            // std::cout << "   " << ilats[I][J] << " " << ilons[I][J] << std::endl;
            if(fabs(ilons[I][J] - currLon) < tol && fabs(ilats[I][J] - currLat) < tol) {
               int startI = std::max(0, I-10);
               int startJ = std::max(0, J-10);
               int endI   = std::min(iFrom.getNumLat()-1, I+10);
               int endJ   = std::min(iFrom.getNumLon()-1, J+10);
               // std::cout << i << " " << j << " " << olats[i][j] << " " << ilats[I][J] << " " << olons[i][j] << " " << ilons[I][J] << std::endl;
               // std::cout << i << " " << j << " " << abs(ilons[I][J] - olons[i][j]) << " " << abs(ilats[I][J] - olats[i][j]) << std::endl;
               // std::cout << i << " " << j << " Searching in [" << startI << " " << startJ << " " << endI << " " << endJ << "]" << std::endl;
               float minDist = Util::MV;
               for(int ii = startI; ii <= endI; ii++) {
                  for(int jj = startJ; jj <= endJ; jj++) {
                     float currDist = Util::getDistance(olats[i][j], olons[i][j], ilats[ii][jj], ilons[ii][jj]);
                     if(!Util::isValid(minDist) || currDist < minDist) {
                        // std::cout << ilats[ii][jj] << " " << ilons[ii][jj] << "    " << currDist << std::endl;
                        I = ii;
                        J = jj;
                        minDist = currDist;
                     }
                  }
               }
               // std::cout << "Found: " << i << " " << j << " " << olats[i][j] << " " << olons[i][j] << " " << ilats[I][J] << " " << ilons[I][J] << " " << minDist << std::endl;
               break;
            }
            else {
               assert(I >= 0);
               assert(J >= 0);
               assert(ilons.size() > I);
               assert(ilons[I].size() > J);
               assert(olons.size() > i);
               assert(olons[i].size() > j);
               assert(ilats.size() > I);
               assert(ilats[I].size() > J);
               assert(olats.size() > i);
               assert(olats[i].size() > j);
               if(ilons[I][J] < currLon-tol)
                  J++;
               else if(ilons[I][J] > currLon+tol)
                  J--;
               if(ilats[I][J] < currLat-tol)
                  I++;
               else if(ilats[I][J] < currLat+tol)
                  I--;
               I = std::min(iFrom.getNumLat()-1, std::max(0, I));
               J = std::min(iFrom.getNumLon()-1, std::max(0, J));
            }
            counter++;
            if(counter > 1000) {
               // std::cout << "Couldn't find for " << i << " " << j << std::endl;
               float minDist = Util::MV;
               for(int ii = 0; ii < iFrom.getNumLat(); ii++) {
                  for(int jj = 0; jj < iFrom.getNumLon(); jj++) {
                     float currDist = Util::getDistance(olats[i][j], olons[i][j], ilats[ii][jj], ilons[ii][jj]);
                     // std::cout << ii << " " << jj << " " << currDist << " " << minDist << std::endl;
                     if(Util::isValid(currDist) && (!Util::isValid(minDist) || currDist < minDist)) {
                        I = ii;
                        J = jj;
                        minDist = currDist;
                     }
                  }
               }
               break;
            }
         }
         iI[i][j] = I;
         iJ[i][j] = J;
         // std::cout << "Found: " << i << " " << j << " " << olats[i][j] << " " << olons[i][j] << " " << ilats[I][J] << " " << ilons[i][J] << std::endl;
      }
   }
   addToCache(iFrom, iTo, iI, iJ);
}

bool Downscaler::isCached(const File& iFrom, const File& iTo) {
   std::map<boost::uuids::uuid, std::map<boost::uuids::uuid, std::pair<vec2Int, vec2Int> > >::const_iterator it = mNeighbourCache.find(iFrom.getUniqueTag());
   if(it == mNeighbourCache.end()) {
      return false;
   }
   std::map<boost::uuids::uuid, std::pair<vec2Int, vec2Int> >::const_iterator it2 = it->second.find(iTo.getUniqueTag());
   if(it2 == it->second.end()) {
      return false;
   }
   return true;
}

void Downscaler::addToCache(const File& iFrom, const File& iTo, vec2Int iI, vec2Int iJ) {
   std::pair<vec2Int, vec2Int> pair(iI, iJ);
   mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()] = pair;
}
bool Downscaler::getFromCache(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(!isCached(iFrom, iTo))
      return false;
   iI = mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()].first;
   iJ = mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()].second;
   return true;
}
