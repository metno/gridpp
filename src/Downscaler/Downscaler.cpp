#include <boost/scoped_ptr.hpp>
#include <cmath>

#include "Downscaler.h"
#include "../File/File.h"

std::map<Uuid, std::map<Uuid, std::pair<vec2Int, vec2Int> > > Downscaler::mNeighbourCache;

namespace {
   class KDTree {
      struct TreeNode {
         bool xsection;
         float lon;
         float lat;
         size_t ipos;
         size_t jpos;

         boost::scoped_ptr<TreeNode> left;
         boost::scoped_ptr<TreeNode> right;

         const TreeNode* parent;

         TreeNode(): xsection(false), left(NULL), right(NULL), parent(NULL) {}
      };

      struct Indexed {
         float lon;
         float lat;
         size_t index1;
         size_t index2;

         Indexed() {}
         Indexed(const float lon_,
                 const float lat_,
                 const size_t index1_,
                 const size_t index2_):
            lon(lon_), lat(lat_), index1(index1_), index2(index2_) {}
      };

      typedef boost::scoped_ptr<TreeNode> unode;
      typedef std::vector<Indexed> indexdVec;

      unode root;

      static const TreeNode* nearestNeighbour(const unode& root, const float lon, const float lat);
      static const TreeNode* firstGuess(const unode& root, const float lon, const float lat);
      static void subTree(indexdVec& iLonLat,
                          const size_t from,
                          const size_t to,
                          const bool xsection,
                          const TreeNode* parent,
                          unode& root);

   public:
      KDTree(): root(NULL) {}
      //KDTree(const KDTree &) = delete;
      //KDTree& operator=(const KDTree &) = delete;

      void buildTree(const File& iFrom);
      void getNearestNeighbour(const File& iTo, vec2Int& iI, vec2Int& iJ) const;

      friend bool compareLons (const KDTree::Indexed& l, const KDTree::Indexed& r);
      friend bool compareLats (const KDTree::Indexed& l, const KDTree::Indexed& r);
   };

   void KDTree::buildTree(const File& iFrom) {

      vec2 ilats = iFrom.getLats();
      vec2 ilons = iFrom.getLons();

      size_t nLon = iFrom.getNumLon();
      size_t nLat = iFrom.getNumLat();

      indexdVec lons(nLon*nLat);

      indexdVec::iterator currLon = lons.begin();
      size_t to = -1;
      for(size_t i = 0; i < ilats.size(); ++i) {
        for(size_t j = 0; j < ilats[0].size(); ++j) {
           if(Util::isValid(ilons[i][j]) && Util::isValid(ilats[i][j])) {
              *(currLon++) = Indexed(ilons[i][j], ilats[i][j], i, j);
              ++to;
           }
        }
      }

      if(to >= 0) subTree(lons, 0, to, true, NULL, root);
   }

   bool compareLons (const KDTree::Indexed& l, const KDTree::Indexed& r){ return (l.lon < r.lon); }
   bool compareLats (const KDTree::Indexed& l, const KDTree::Indexed& r){ return (l.lat < r.lat); }

   void KDTree::subTree(indexdVec& iLonLat,
                        const size_t from,
                        const size_t to,
                        const bool xsection,
                        const TreeNode* parent,
                        unode& root) {

      unode node(new KDTree::TreeNode());
      root.swap(node);

      root->parent = parent;
      root->xsection = xsection;

      size_t len = to - from + 1;

      if(len == 1) {
         root->lon = iLonLat[from].lon;
         root->lat = iLonLat[from].lat;
         root->ipos = iLonLat[from].index1;
         root->jpos = iLonLat[from].index2;
         return;
      }

      size_t med = from + len/2;
      if(xsection) {
         std::sort(iLonLat.begin() + from, iLonLat.begin() + to + 1, compareLons);
         while(iLonLat[med].lon == iLonLat[med+1].lon && med < to) { ++med; }
      } else {
         std::sort(iLonLat.begin() + from, iLonLat.begin() + to + 1, compareLats);
         while(iLonLat[med].lat == iLonLat[med+1].lat && med < to) { ++med; }
      }

      root->lon = iLonLat[med].lon;
      root->lat = iLonLat[med].lat;
      root->ipos = iLonLat[med].index1;
      root->jpos = iLonLat[med].index2;

      if(med > from) subTree(iLonLat, from, med - 1, !xsection, root.get(), root->left);
      if(med < to) subTree(iLonLat, med + 1, to, !xsection, root.get(), root->right);

   }

   const KDTree::TreeNode* KDTree::firstGuess(const unode& root, const float lon, const float lat) {
      if(root->xsection){
         if(lon <= root->lon) {
            if(root->left) return firstGuess(root->left, lon, lat);
         } else {
            if(root->right) return firstGuess(root->right, lon, lat);
         }
      } else {
         if(lat <= root->lat) {
            if(root->left) return firstGuess(root->left, lon, lat);
         } else {
            if(root->right) return firstGuess(root->right, lon, lat);
         }
      }
      return root.get();
   }

   const KDTree::TreeNode* KDTree::nearestNeighbour(const unode& root, const float lon, const float lat) {

      const TreeNode * nearLeaf = firstGuess(root, lon, lat);

      float currBest = Util::getDistance(lat, lon, nearLeaf->lat, nearLeaf->lon);

      const TreeNode * currentLeaf = nearLeaf;
      const TreeNode * nearestLeaf = nearLeaf;

      while (currentLeaf) {
         float dist = Util::getDistance(lat, lon, currentLeaf->lat, currentLeaf->lon);

         if(dist < currBest) {
            currBest = dist;
            nearestLeaf = currentLeaf;
         }

         bool cross;
         if(currentLeaf->xsection) {
            cross = Util::getDistance(lat, lon, lat, currentLeaf->lon) < currBest;
         } else {
            cross = Util::getDistance(lat, lon, currentLeaf->lat, lon) < currBest;
         }

         if(cross) {
            if(currentLeaf->left.get() == nearLeaf) {
               nearLeaf = NULL;
               if(currentLeaf->right) nearLeaf = nearestNeighbour( currentLeaf->right, lon, lat);
            } else {
               nearLeaf = NULL;
               if(currentLeaf->left) nearLeaf = nearestNeighbour(currentLeaf->left, lon, lat);
            }

            if(nearLeaf) {
               float dist = Util::getDistance(lat, lon, nearLeaf->lat, nearLeaf->lon);
               if(dist < currBest) {
                  currBest = dist;
                  nearestLeaf = nearLeaf;
               }
            }
         }

         nearLeaf = currentLeaf;
         if(currentLeaf == root.get()) {
            currentLeaf = NULL;
         } else {
            currentLeaf = currentLeaf->parent;
         }
      }

      return nearestLeaf;
   }

   void KDTree::getNearestNeighbour(const File& iTo, vec2Int& iI, vec2Int& iJ) const {

      vec2 olats = iTo.getLats();
      vec2 olons = iTo.getLons();
      size_t nLon = iTo.getNumLon();
      size_t nLat = iTo.getNumLat();

      iI.resize(nLat);
      iJ.resize(nLat);

      if(!root) {
         for(size_t i = 0; i < nLat; ++i) {
            iI[i].resize(nLon, Util::MV);
            iJ[i].resize(nLon, Util::MV);
         }
         return;
      }

      const TreeNode * nearest;

      #pragma omp parallel for private(nearest)
      for(size_t i = 0; i < nLat; ++i) {
         iI[i].resize(nLon, Util::MV);
         iJ[i].resize(nLon, Util::MV);
         for(size_t j = 0; j < nLon; ++j) {
            if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
               // Find the nearest neighbour from input grid (ii, jj)
                  nearest = nearestNeighbour(root, olons[i][j], olats[i][j]);
                  iI[i][j] = nearest->ipos;
                  iJ[i][j] = nearest->jpos;
            }
         }
      }
   }

}

Downscaler::Downscaler(Variable::Type iVariable, const Options& iOptions) : Scheme(iOptions),
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
      return new DownscalerNearestNeighbour(iVariable, iOptions);
   }
   else if(iName == "gradient") {
      DownscalerGradient* d = new DownscalerGradient(iVariable, iOptions);
      return d;
   }
   else if(iName == "smart") {
      DownscalerSmart* d = new DownscalerSmart(iVariable, iOptions);
      float searchRadius = Util::MV;
      if(iOptions.getValue("searchRadius", searchRadius)) {
         d->setSearchRadius(searchRadius);
      }
      float numSmart = Util::MV;
      if(iOptions.getValue("numSmart", numSmart)) {
         d->setNumSmart(numSmart);
      }
      float minElevDiff = Util::MV;
      if(iOptions.getValue("minElevDiff", minElevDiff)) {
         d->setMinElevDiff(minElevDiff);
      }
      return d;
   }
   else if(iName == "bypass") {
      DownscalerBypass* d = new DownscalerBypass(iVariable, iOptions);
      return d;
   }
   else if(iName == "pressure") {
      DownscalerPressure* d = new DownscalerPressure(iVariable, iOptions);
      return d;
   }
   else {
      Util::error("Could not instantiate downscaler of type '" + iName + "'");
      return NULL;
   }
}

void Downscaler::getNearestNeighbourBruteForce(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
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
         Util::status("Grids are identical, short cut in finding nearest neighbours");
         addToCache(iFrom, iTo, iI, iJ);
         return;
      }
   }

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      iI[i].resize(nLon, Util::MV);
      iJ[i].resize(nLon, Util::MV);
      for(int j = 0; j < nLon; j++) {
         if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
            getNearestNeighbourBruteForce(iFrom, olons[i][j], olats[i][j], iI[i][j], iJ[i][j]);
         }
      }
   }
   addToCache(iFrom, iTo, iI, iJ);
}

void Downscaler::getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(isCached(iFrom, iTo)) {
       getFromCache(iFrom, iTo, iI, iJ);
       return;
   }

   // Check if the grid is the same
   if(iFrom.getNumLat() == iTo.getNumLat() && iFrom.getNumLon() == iTo.getNumLon()) {
      vec2 ilats = iFrom.getLats();
      vec2 ilons = iFrom.getLons();
      vec2 olats = iTo.getLats();
      vec2 olons = iTo.getLons();

      if(ilats == olats && ilons == olons) {
         int nLon = iTo.getNumLon();
         int nLat = iTo.getNumLat();

         iI.resize(nLat);
         iJ.resize(nLat);

         for(int i = 0; i < nLat; i++) {
            iI[i].resize(nLon, 0);
            iJ[i].resize(nLon, 0);
            for(int j = 0; j < nLon; j++) {
               iI[i][j] = i;
               iJ[i][j] = j;
            }
         }
         Util::status("Grids are identical, short cut in finding nearest neighbours");
         addToCache(iFrom, iTo, iI, iJ);
         return;
      }
   }

   KDTree searchTree;
   searchTree.buildTree(iFrom);
   searchTree.getNearestNeighbour(iTo, iI, iJ);

   addToCache(iFrom, iTo, iI, iJ);
}

bool Downscaler::isCached(const File& iFrom, const File& iTo) {
   std::map<Uuid, std::map<Uuid, std::pair<vec2Int, vec2Int> > >::const_iterator it = mNeighbourCache.find(iFrom.getUniqueTag());
   if(it == mNeighbourCache.end()) {
      return false;
   }
   std::map<Uuid, std::pair<vec2Int, vec2Int> >::const_iterator it2 = it->second.find(iTo.getUniqueTag());
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

void Downscaler::getNearestNeighbourBruteForce(const File& iFrom, float iLon, float iLat, int& iI, int &iJ) {
   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();

   iI = Util::MV;
   iJ = Util::MV;

   if(Util::isValid(iLat) && Util::isValid(iLon)) {
      float minDist = Util::MV;
      for(int i = 0; i < ilats.size(); i++) {
         for(int j = 0; j < ilats[0].size(); j++) {
            if(Util::isValid(ilats[i][j]) && Util::isValid(ilons[i][j])) {
               float currDist = Util::getDistance(ilats[i][j], ilons[i][j], iLat, iLon);
               if(!Util::isValid(minDist) || currDist < minDist) {
                  iI = i;
                  iJ = j;
                  minDist = currDist;
               }
            }
         }
      }
   }
}

std::string Downscaler::getDescriptions() {
   std::stringstream ss;
   ss << DownscalerNearestNeighbour::description();
   ss << DownscalerGradient::description();
   ss << DownscalerSmart::description();
   ss << DownscalerPressure::description();
   ss << DownscalerBypass::description();
   return ss.str();
}
