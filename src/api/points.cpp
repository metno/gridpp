#include "gridpp.h"

gridpp::Points::Points() {
    // TODO: Deal with the empty case. Can't really find nearest neighbours then

}
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
ivec gridpp::Points::get_in_domain_indices(const gridpp::Grid& grid) const {
    ivec indices;
    indices.reserve(size());
    ivec gsize = grid.size();
    int Y = gsize[0];
    int X = gsize[1];

    for(int s = 0; s < size(); s++) {
        const std::vector<int> ind = grid.get_nearest_neighbour(mLats[s], mLons[s]);
        int y = ind[0];
        int x = ind[1];
        if(x > 0 && x < X-1 && y > 0 && y < Y-1) {
            indices.push_back(s);
        }
    }
    return indices;
}
gridpp::Points gridpp::Points::get_in_domain(const gridpp::Grid& grid) const {
    vec lats, lons, elevs, lafs;
    lats.reserve(size());
    lons.reserve(size());
    elevs.reserve(size());
    lafs.reserve(size());
    ivec indices = get_in_domain_indices(grid);
    int S = indices.size();
    for(int s = 0; s < S; s++) {
        lats.push_back(mLats[indices[s]]);
        lons.push_back(mLons[indices[s]]);
        elevs.push_back(mElevs[indices[s]]);
        lafs.push_back(mLafs[indices[s]]);
    }
    const gridpp::Points new_points(lats, lons, elevs, lafs);
    return new_points;
}
gridpp::Points& gridpp::Points::operator=(gridpp::Points other) {
    std::swap(mLats, other.mLats);
    std::swap(mLons, other.mLons);
    std::swap(mElevs, other.mElevs);
    std::swap(mLafs, other.mLafs);
    std::swap(mTree, other.mTree);
    return *this;
}
gridpp::Points::Points(const gridpp::Points& other) {
    mLats = other.mLats;
    mLons = other.mLons;
    mElevs = other.mElevs;
    mLafs = other.mLafs;
    mTree = KDTree(mLats, mLons);
}
