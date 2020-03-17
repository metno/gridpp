#include "gridpp.h"

gridpp::Points::Points(vec lats, vec lons, vec elevs, vec lafs) {
    mLats = lats;
    mLons = lons;
    mElevs = elevs;
    mLafs = lafs;
    KDTree tree = KDTree(lats, lons);
    mTree = tree;
}

int gridpp::Points::get_num_neighbours(float lat, float lon, float radius) {
    return mTree.get_num_neighbours(lat, lon, radius);
}

ivec gridpp::Points::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) {
    return mTree.get_neighbours_with_distance(lat, lon, radius, distances);
}

ivec gridpp::Points::get_neighbours(float lat, float lon, float radius) {
    return mTree.get_neighbours(lat, lon, radius);
}

ivec gridpp::Points::get_closest_neighbours(float lat, float lon, int num) {
    return mTree.get_closest_neighbours(lat, lon, num);
}
int gridpp::Points::get_nearest_neighbour(float lat, float lon) {
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
