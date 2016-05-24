#include "VPTree.h"
#include <algorithm>
#include <limits>

namespace {
   bool operator==(const Sincos &l, const Sincos &r) {
      return l.cosine == r.cosine && r.sine == l.sine;
   }
}

void VPTree::subTree(const size_t from, const size_t to, VPTree::unode& root) {

   unode node(new VPTree::TreeNode());
   root.swap(node);

   size_t len = to - from + 1;
   root->index = from;

   if ( len > 1 ) {
      int i = static_cast<int>(static_cast<double>(rand()) / RAND_MAX * (to - from - 1)) + from;

      std::swap( mCoords[from], mCoords[i] );

      int median = from + len/2;

      double medianNextDistance = 0;
      if(len == 2){
         median = from;
         medianNextDistance = VPTree::getDistance(
                     mCoords[from].lat, mCoords[from].lon,
                     mCoords[from+1].lat, mCoords[from+1].lon );
      } else {
         std::nth_element(
                     mCoords.begin() + from + 1,
                     mCoords.begin() + median,
                     mCoords.begin() + to + 1,
                     DistanceComparator( mCoords[from] ));
         if(median%2 == 0) {
             std::vector<Indexed>::iterator medianNext = std::min_element(
                         mCoords.begin() + median + 1,
                         mCoords.begin() + to + 1,
                         DistanceComparator( mCoords[from] ));
             medianNextDistance = VPTree::getDistance(
                         mCoords[from].lat, mCoords[from].lon,
                         medianNext->lat, medianNext->lon );
         }
      }

      double medianDistance = VPTree::getDistance(
                 mCoords[from].lat, mCoords[from].lon,
                 mCoords[median].lat, mCoords[median].lon );

      root->cutDistance = medianDistance + medianNextDistance;
      if(median%2 == 0) root->cutDistance /= 2.;


#pragma omp task shared(root)
      subTree( from + 1, median, root->left );
#pragma omp task shared(root)
      if(median < to) subTree( median + 1, to, root->right );
#pragma omp taskwait
   }
}

void VPTree::nearestNeighbour(const VPTree::unode& root, const double lon, const double lat, double& minDist, size_t& bestId) const
{
   double currDist = VPTree::getDistance(lat, lon, mCoords[root->index].lat, mCoords[root->index].lon);

   if(currDist < minDist) {
      minDist = currDist;
      bestId = root->index;
   }

   if(currDist < root->cutDistance) {
      if (root->left) {
         nearestNeighbour(root->left, lon, lat, minDist, bestId);
      }

      if (currDist + minDist >= root->cutDistance && root-> right) {
         nearestNeighbour( root->right, lon, lat, minDist, bestId );
      }
   } else {
      if (root->right) {
         nearestNeighbour(root->right, lon, lat, minDist, bestId);
      }

      if (currDist - minDist <= root->cutDistance && root->left) {
         nearestNeighbour(root->left, lon, lat, minDist, bestId);
      }
   }

}

void VPTree::build(const vec2& iLats, const vec2& iLons)
{

   if(iLats.size() != iLons.size())
      Util::error("Cannot initialize VPTree, lats and lons not the same size");

   size_t nLat = iLats.size();
   if(nLat == 0)
      Util::error("Cannot initialize VPTree, no valid locations");

   size_t nLon = iLats[0].size();
   if(nLon == 0)
      Util::error("Cannot initialize VPTree, no valid locations");

   mCoords.reserve(nLon*nLat);

   size_t to = -1;
   for(size_t i = 0; i < nLat; ++i) {
      for(size_t j = 0; j < nLon; ++j) {
         if(Util::isValid(iLons[i][j]) && Util::isValid(iLats[i][j])) {

            assert(fabs(iLats[i][j]) <= 90 && fabs(iLons[i][j]) <= 360);

            mCoords.push_back(Indexed(iLons[i][j], iLats[i][j], i, j));
            ++to;
         }
      }
   }

   if(to == -1) {
      Util::error("Cannot initialize VPTree, no valid locations");
   }
#pragma omp parallel
#pragma omp single
   {
      if(to >= 0) subTree(0, to, mRoot);
   }
}

void VPTree::getNearestNeighbour(const File& iTo, vec2Int& iI, vec2Int& iJ) const {
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   size_t nLon = iTo.getNumLon();
   size_t nLat = iTo.getNumLat();

   iI.resize(nLat);
   iJ.resize(nLat);

   if(!mRoot) {
      for(size_t i = 0; i < nLat; ++i) {
         iI[i].resize(nLon, Util::MV);
         iJ[i].resize(nLon, Util::MV);
      }
      return;
   }

   const size_t bigVal = mCoords.size();
#pragma omp parallel for
   for(size_t i = 0; i < nLat; ++i) {
      iI[i].resize(nLon, Util::MV);
      iJ[i].resize(nLon, Util::MV);
      for(size_t j = 0; j < nLon; ++j) {
         if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
            // Find the nearest neighbour from input grid (ii, jj)
            double nrstDistance = 1.E20;
            size_t nrstIndex = bigVal;

            nearestNeighbour(mRoot, olons[i][j], olats[i][j], nrstDistance, nrstIndex);
            if(nrstIndex < bigVal) {
               iI[i][j] = mCoords[nrstIndex].ipos;
               iJ[i][j] = mCoords[nrstIndex].jpos;
            }
         }
      }
   }
}

double VPTree::getDistance(const Sincos& lat1, const Sincos& lon1, const Sincos& lat2, const Sincos& lon2) {
   if(lat1 == lat2 && lon1 == lon2)
      return 0;

   double ratio = lat1.cosine*lon1.cosine*lat2.cosine*lon2.cosine
                + lat1.cosine*lon1.sine*lat2.cosine*lon2.sine
                + lat1.sine*lat2.sine;
   double dist = acos(ratio)*Util::radiusEarth;

   return dist;
}
