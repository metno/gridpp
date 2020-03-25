#include "gridpp.h"

gridpp::Points::Points(vec lats, vec lons, vec elevs, vec lafs) {
    mLats = lats;
    mLons = lons;
    mElevs = elevs;
    mLafs = lafs;
    KDTree tree = KDTree(lats, lons);
    mTree = tree;
    int N = mLats.size();
    if(mElevs.size() != N) {
        mElevs.clear();
        for(int i = 0; i < N; i++)
            mElevs.push_back(0);
    }
    if(mLafs.size() != N) {
        mLafs.clear();
        for(int i = 0; i < N; i++)
            mLafs.push_back(1);
    }
}

int gridpp::Points::get_num_neighbours(float lat, float lon, float radius) const {
    return mTree.get_num_neighbours(lat, lon, radius);
}

ivec gridpp::Points::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const {
    return mTree.get_neighbours_with_distance(lat, lon, radius, distances);
}

ivec gridpp::Points::get_neighbours(float lat, float lon, float radius) const {
    return mTree.get_neighbours(lat, lon, radius);
}

ivec gridpp::Points::get_closest_neighbours(float lat, float lon, int num) const {
    return mTree.get_closest_neighbours(lat, lon, num);
}
int gridpp::Points::get_nearest_neighbour(float lat, float lon) const {
    return get_closest_neighbours(lat, lon, 1)[0];
}
vec gridpp::Points::get_lats() const {
    return mLats;
}
vec gridpp::Points::get_lons() const {
    return mLons;
}
vec gridpp::Points::get_elevs() const {
    return mElevs;
}
vec gridpp::Points::get_lafs() const {
    return mLafs;
}
int gridpp::Points::size() const {
    return mLats.size();
}
