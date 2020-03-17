#ifndef API_H
#define API_H
#include <vector>
#include <set>
#include <string>
#include <armadillo>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>
typedef std::vector<std::vector<std::vector<float> > > vec3;
typedef std::vector<std::vector<float> > vec2;
typedef std::vector<float> vec;
// typedef std::vector<float> fvec;
typedef std::vector<int> ivec;
namespace gridpp {
    typedef arma::mat mattype;
    typedef arma::vec vectype;
    typedef arma::cx_mat cxtype;

    float getDistance(float lat1, float lon1, float lat2, float lon2, bool approx=false);
    bool isValid(float iValue);
    float deg2rad(float deg);
    float rad2deg(float rad);
    void debug(std::string string);
    void error(std::string string);
    float calcRho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wmin);

    static float MV;
    static float pi;
    static double radiusEarth;

    void check_vec(vec2 input, int Y, int X);
    void check_vec(vec input, int S);

    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };

    class KDTree {
        public:
            KDTree(vec lats, vec lons, vec elevs=vec(), vec lafs=vec());
            KDTree(vec2 lats, vec2 lons, vec2 elevs=vec2(), vec2 lafs=vec2());
            KDTree(const KDTree& other) {abort();};
            KDTree& operator=(KDTree other) {abort();};

            /** Find single nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             * */
            int get_nearest_neighbour(float lat, float lon);
            ivec get_nearest_neighbour_2d(float lat, float lon);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            ivec get_neighbours(float lat, float lon, float radius);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             * */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances);

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            int get_num_neighbours(float lat, float lon, float radius);

            /** Find a set of nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param num Number of points to find
             * */
            ivec get_closest_neighbours(float lat, float lon, int num);


            /** Convert lat/lons to 3D cartesian coordinates with the centre of the earth as the origin
             *  @param lats vector of latitudes [deg]
             *  @param lons vector of longitudes [deg]
             *  @param x_coords vector of x-coordinates [m]
             *  @param y_coords vector of y-coordinates [m]
             *  @param z_coords vector of z-coordinates [m]
             * */
            static bool convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords);

            /** Same as above, but convert a single lat/lon to 3D cartesian coordinates
             *  @param lat latitude [deg]
             *  @param lon longitude [deg]
             *  @param x_coord x-coordinate [m]
             *  @param y_coord y-coordinate [m]
             *  @param z_coord z-coordinate [m]
             * */
            static bool convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord);
            static float deg2rad(float deg);
            static float calc_distance(float lat1, float lon1, float lat2, float lon2);
            static float calc_distance(float x0, float y0, float z0, float x1, float y1, float z1);
            vec get_lats() const;
            vec get_lons() const;
            vec get_elevs() const;
            vec get_lafs() const;
            vec2 get_lats_2d() const;
            vec2 get_lons_2d() const;
            vec2 get_elevs_2d() const;
            vec2 get_lafs_2d() const;
        private:
            typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
            typedef std::pair<point, unsigned> value;
            typedef boost::geometry::model::box<point> box;
            boost::geometry::index::rtree< value, boost::geometry::index::quadratic<16> > mTree;
            int mX;
            vec mLats;
            vec mLons;
            vec mElevs;
            vec mLafs;
            vec2 get_2d(vec input) const;
    };

    int optimal_interpolation_single_member(const vec2& input,
            const gridpp::KDTree& bgrid,
            const vec& pobs,
            const vec& pci,
            const gridpp::KDTree& ptree,
            float minRho,
            float hlength,
            float vlength,
            float wmin,
            float maxElevDiff,
            bool landOnly,
            int maxLocations,
            float elevGradient,
            float epsilon,
            vec2& output);

    int optimal_interpolation_single_member(const vec2& input,
            const vec2& blats,
            const vec2& blons,
            const vec2& belevs,
            const vec2& blafs,
            const vec& pobs,
            const vec& pci,
            const vec& plats,
            const vec& plons,
            const vec& pelevs,
            const vec& plafs,
            float minRho,
            float hlength,
            float vlength,
            float wmin,
            float maxElevDiff,
            bool landOnly,
            int maxLocations,
            float elevGradient,
            float epsilon,
            vec2& output);

    vec2 neighbourhood(const vec2& input, int radius, std::string operation, float quantile=-1);
    vec3 neighbourhood(const vec3& input, int radius, std::string operation, float quantile=-1);
    // vec2 bilinear(const KDTree& tree, const vec2

    namespace util {
        bool is_valid(float value);
        float calculate_stat(const std::vector<float>& iArray, std::string iOperator, float iQuantile=MV);
        int num_missing_values(const vec2& iArray);
    }
    /*
    namespace downscaling {
        vec2 nearestNeighbour(const vec2& iLats, const vec2& iLons, const vec2& iValues, const vec2& oLats, const vec2& oLons, int iRadius, std::string iStatType, std::string iQuantile);
        vec2 nearestNeighbour(Field ifield, Grid ogrid, int iRadius, std::string iStatType, std::string iQuantile);
    }
    namespace calibration {
        vec2 quantileMapping(const vec2& iValues, std::string iExtrapolation, vec iQuantiles, const ParameterSet& iParameterSet);
        vec2 neighbourhood(const Field& iField, const Grid& oGrid, float iRadius, std::string iOperator);
        vec2 optimal_interpolation_single_member(const vec2& iValues, const Grid& iGrid);
        vec2 optimal_interpolation_ensemble(const vec3& iValues, const Grid& iGrid, const ParameterSet& iParameterSet);
    }
    namespace diagnose {
        vec2 computeRh(const vec2& iTemperatures, const vec2& iDewPointTemperatures);
        float computeRh(float iTemperatures, float iDewPointTemperatures);
    }

    class Grid(Dataset) {
        Grid(vec2 lats, vec2 lons, vec2 elevs, vec2 lafs);
    }
    class ParameterSet {
        ParameterSetSpatial(vec2 iLats, vec2 iLons, vec3 iRef, vec3 iModel);
        void get(float iLat, float iLon, float iElev, vec iRef, vec iModel);
        private:
            vec2 iLats, iLons;
            vec3 iRef, iModel;
    };
    */
};
#endif
