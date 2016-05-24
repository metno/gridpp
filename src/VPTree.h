#ifndef VPTREE_H
#define VPTREE_H
#include <boost/scoped_ptr.hpp>
#include <cmath>
#include <stack>
#include <vector>
#include "Util.h"
#include "File/File.h"

typedef std::vector<std::vector<int> > vec2Int;

struct Sincos {
   double sine;
   double cosine;
   Sincos(double val) {
      sine = sin(Util::deg2rad(val));
      cosine = cos(Util::deg2rad(val));
   }
};

class VPTree {
   struct TreeNode {
      size_t index;
      double cutDistance;

      boost::scoped_ptr<TreeNode> left;
      boost::scoped_ptr<TreeNode> right;

      TreeNode(): index(0), cutDistance(0.),left(NULL), right(NULL) {}
   };

   struct Indexed {
      Sincos lon;
      Sincos lat;
      size_t ipos;
      size_t jpos;

      Indexed(const double lon_,
              const double lat_,
              const size_t index1_,
              const size_t index2_):
         lon(lon_), lat(lat_), ipos(index1_), jpos(index2_) {}
   };

   struct DistanceComparator {
      const Indexed& base;
      DistanceComparator( const Indexed& base_ ) : base(base_) {}
      bool operator()(const Indexed& l, const Indexed& r) {
         double dl = VPTree::getDistance(base.lat, base.lon, l.lat, l.lon);
         double dr = VPTree::getDistance(base.lat, base.lon, r.lat, r.lon);
         return dl < dr;
      }
    };

   typedef boost::scoped_ptr<TreeNode> unode;

   std::vector<Indexed> mCoords;
   unode mRoot;

   void subTree(const size_t from, const size_t to, unode& root);

   void nearestNeighbour(const unode& root, const double lon, const double lat, double& minDist, size_t &bestId ) const;

public:
   VPTree():mRoot(NULL){}
   VPTree(const vec2& iLats, const vec2& iLons){
      build(iLats, iLons);
   }

   void build(const vec2& iLats, const vec2& iLons);

   void getNearestNeighbour(const File& iTo, vec2Int& iI, vec2Int& iJ) const;

   ~VPTree() {}

   static double getDistance(const Sincos& lat1, const Sincos& lon1, const Sincos& lat2, const Sincos& lon2);
};

#endif // VPTREE_H
