#include "gridpp.h"
#include <iostream>

using namespace gridpp;

gridpp::KDTree::KDTree(vec lats, vec lons, CoordinateType type) {
    mLats = lats;
    mLons = lons;
    mType = type;

    gridpp::convert_coordinates(mLats, mLons, mX, mY, mZ, mType);
    for(int i = 0; i < mLats.size(); i++) {
        point p(mX[i], mY[i], mZ[i]);
        mTree.insert(std::make_pair(p, i));
    }
}

int gridpp::KDTree::get_num_neighbours(float lat, float lon, float radius, bool include_match) const {
    ivec indices = get_neighbours(lat, lon, radius, include_match);
    return indices.size();
}

ivec gridpp::KDTree::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match) const {
    float x, y, z;
    gridpp::convert_coordinates(lat, lon, x, y, z, mType);
    ivec indices = get_neighbours(lat, lon, radius, include_match);

    int num = indices.size();
    distances.resize(num);
    for(int i = 0; i < num; i++) {
        float x1, y1, z1;
        gridpp::convert_coordinates(mLats[indices[i]], mLons[indices[i]], x1, y1, z1, mType);
        distances[i] = gridpp::KDTree::calc_straight_distance(x, y, z, x1, y1, z1);
    }

    return indices;
}

ivec gridpp::KDTree::get_neighbours(float lat, float lon, float radius, bool include_match) const {
    float x, y, z;
    gridpp::convert_coordinates(lat, lon, x, y, z, mType);

    std::vector<value> results;
#if 1
    point p(x, y, z);
    box bx(point(x - radius, y - radius, z - radius), point(x + radius, y + radius, z + radius));
    within_radius r(p, radius, include_match);

    // NOTE: The boost::geometry::index::within predicate is extremely efficient at finding points
    // within the box, whereas the within_radius predicate is slow since boost has to test all points.
    // By including within_radius second in the && only when the within predicate evaluates to true
    // will the within_radius predicate be evaluated. See #6e4d5ca for how to do this incorrectly.
    mTree.query(boost::geometry::index::within(bx) && boost::geometry::index::satisfies(r), std::back_inserter(results));
    int num = results.size();

    ivec ret;
    ret.reserve(num);
    for(int i = 0; i < num; i++) {
        ret.push_back(results[i].second);
    }
#else
    // Alternative implementation that is roughly the same speed
    box bx(point(x - radius, y - radius, z - radius), point(x + radius, y + radius, z + radius));
    mTree.query(boost::geometry::index::within(bx), std::back_inserter(results));
    int num = results.size();

    ivec ret;
    ret.reserve(num);
    for(int i = 0; i < num; i++) {
        float x1, y1, z1;
        gridpp::convert_coordinates(mLats[results[i].second], mLons[results[i].second], x1, y1, z1, mType);
        float dist = gridpp::KDTree::calc_straight_distance(x, y, z, x1, y1, z1);
        if(dist <= radius) {
            if(include_match || dist != 0)
                ret.push_back(results[i].second);
        }
    }
#endif
    return ret;
}

ivec gridpp::KDTree::get_closest_neighbours(float lat, float lon, int num, bool include_match) const {
    float x, y, z;
    gridpp::convert_coordinates(lat, lon, x, y, z, mType);
    point p(x, y, z);

    std::vector<value> results;
    if(!include_match) {
        is_not_equal s(p);
        mTree.query(boost::geometry::index::nearest(p, num) && boost::geometry::index::satisfies(s), std::back_inserter(results));
    }
    else {
        mTree.query(boost::geometry::index::nearest(p, num), std::back_inserter(results));
    }
    int num_found = results.size();

    ivec ret;
    ret.reserve(num);
    for(int i = 0; i < num_found; i++) {
        ret.push_back(results[i].second);
    }
    return ret;
}
int gridpp::KDTree::get_nearest_neighbour(float lat, float lon, bool include_match) const {
    return get_closest_neighbours(lat, lon, 1, include_match)[0];
}
float gridpp::KDTree::calc_distance(float lat1, float lon1, float lat2, float lon2, CoordinateType type) {
    if(type == gridpp::Cartesian) {
        float dx = lon1 - lon2;
        float dy = lat1 - lat2;
        return sqrt(dx * dx + dy * dy);
    }
    else if(type == gridpp::Geodetic){
        if(!(fabs(lat1) <= 90 && fabs(lat2) <= 90 && fabs(lon1) <= 360 && fabs(lon2) <= 360)) {
            // std::cout << " Cannot calculate distance, invalid lat/lon: (" << lat1 << "," << lon1 << ") (" << lat2 << "," << lon2 << ")";
            // std::cout << '\n';
        }
        if(lat1 == lat2 && lon1 == lon2)
            return 0;

        double lat1r = deg2rad(lat1);
        double lat2r = deg2rad(lat2);
        double lon1r = deg2rad(lon1);
        double lon2r = deg2rad(lon2);
        double radiusEarth = 6.378137e6;

        double ratio = cos(lat1r)*cos(lon1r)*cos(lat2r)*cos(lon2r)
                       + cos(lat1r)*sin(lon1r)*cos(lat2r)*sin(lon2r)
                       + sin(lat1r)*sin(lat2r);
        double dist = acos(ratio)*radiusEarth;
        return (float) dist;
    }
}
float gridpp::KDTree::calc_distance_fast(float lat1, float lon1, float lat2, float lon2, CoordinateType type) {
    if(type == gridpp::Cartesian) {
        float dx = lon1 - lon2;
        float dy = lat1 - lat2;
        return sqrt(dx * dx + dy * dy);
    }
    else if(type == gridpp::Geodetic){
        double lat1r = deg2rad(lat1);
        double lat2r = deg2rad(lat2);
        double lon1r = deg2rad(lon1);
        double lon2r = deg2rad(lon2);

        /*
        // Convert to 3D coordinates and calculate distance
        // std::cout << lon << " " << lat << std::endl;
        double x1 = std::cos(lat1r) * std::cos(lon1r);
        double y1 = std::cos(lat1r) * std::sin(lon1r);
        float z1 = std::sin(lat1r);

        double x2 = std::cos(lat2r) * std::cos(lon2r);
        double y2 = std::cos(lat2r) * std::sin(lon2r);
        float z2 = std::sin(lat2r);

        float dx = x1 - x2;
        float dy = y1 - y2;
        float dz = z1 - z2;
        return sqrt(dx*dx + dy*dy + dz*dz) * gridpp::radius_earth;
        */

        double dlon = fmod(fabs(lon1r - lon2r), 2 * M_PI);
        // Deal with the wrapping of the longitudes
        if(dlon > M_PI)
            dlon = 2 * M_PI - dlon;

        double mean_lat = (lat1r + lat2r) / 2;
        // Using the mean latitude to account for narrowing of longitude separation doesn't work
        // so well near the poles. Try using the maximum absolute vaulue of the latitude.
        double max_lat = lat1r;
        if (fabs(lat2r) > fabs(lat1r))
            max_lat = lat2r;

        float dx2 = pow(cos(max_lat), 2) * dlon * dlon;
        float dy2 = (lat1r - lat2r) * (lat1r - lat2r);
        return gridpp::radius_earth * sqrt(dx2 + dy2);
    }
    else {
        throw std::runtime_error("Unknown coordinate type");
    }
}
float gridpp::KDTree::calc_distance(const Point3D& p1, const Point3D& p2) {
    if(p1.type != p2.type)
        throw std::runtime_error("Coordinate types must be the same");

    return calc_distance(p1.lat, p1.lon, p2.lat, p2.lon, p1.type);
}
float gridpp::KDTree::calc_straight_distance(const Point3D& p1, const Point3D& p2) {
    return calc_straight_distance(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
}
float gridpp::KDTree::calc_straight_distance(float x0, float y0, float z0, float x1, float y1, float z1) {
    return sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1));
}
float gridpp::KDTree::deg2rad(float deg) {
   return (deg * M_PI / 180);
}
float gridpp::KDTree::rad2deg(float rad) {
   return (rad * 180 / M_PI);
}
vec gridpp::KDTree::get_lats() const {
    return mLats;
}
vec gridpp::KDTree::get_lons() const {
    return mLons;
}
const vec& gridpp::KDTree::get_x() const {
    return mX;
}
const vec& gridpp::KDTree::get_y() const {
    return mY;
}
const vec& gridpp::KDTree::get_z() const {
    return mZ;
}
int gridpp::KDTree::size() const {
    return mLats.size();
}
CoordinateType gridpp::KDTree::get_coordinate_type() const {
    return mType;
}
gridpp::KDTree& gridpp::KDTree::operator=(gridpp::KDTree other) {
    std::swap(mX, other.mX);
    std::swap(mY, other.mY);
    std::swap(mZ, other.mZ);
    std::swap(mLats, other.mLats);
    std::swap(mLons, other.mLons);
    std::swap(mTree, other.mTree);
    std::swap(mType, other.mType);
    return *this;
}
gridpp::KDTree::KDTree(const gridpp::KDTree& other) {
    mX = other.mX;
    mY = other.mY;
    mZ = other.mZ;
    mLats = other.mLats;
    mLons = other.mLons;
    mTree = other.mTree;
    mType = other.mType;
}
gridpp::KDTree::within_radius::within_radius(point p, float radius, bool include_match)  {
    this->p = p;
    this->radius = radius;
    this->include_match = include_match;
}

bool gridpp::KDTree::within_radius::operator()(value const& v) const {
    float x0 = v.first.get<0>();
    float y0 = v.first.get<1>();
    float z0 = v.first.get<2>();
    float x1 = p.get<0>();
    float y1 = p.get<1>();
    float z1 = p.get<2>();
    if(include_match)
        return gridpp::KDTree::calc_straight_distance(x0, y0, z0, x1, y1, z1) <= radius;
    else {
        float dist = gridpp::KDTree::calc_straight_distance(x0, y0, z0, x1, y1, z1);
        return dist <= radius && dist > 0;
    }
}
gridpp::KDTree::is_not_equal::is_not_equal(point p)  {
    this->p = p;
}

bool gridpp::KDTree::is_not_equal::operator()(value const& v) const {
    float x0 = v.first.get<0>();
    float y0 = v.first.get<1>();
    float z0 = v.first.get<2>();
    return p.get<0>() != x0 || p.get<1>() != y0 || p.get<2>() != z0;
}
