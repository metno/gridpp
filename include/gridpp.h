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
#define GRIDPP_VERSION "0.3.2-dev"
#define __version__ GRIDPP_VERSION
typedef std::vector<std::vector<std::vector<float> > > vec3;
typedef std::vector<std::vector<float> > vec2;
typedef std::vector<float> vec;
// typedef std::vector<float> fvec;
typedef std::vector<int> ivec;
typedef std::vector<std::vector<int> > ivec2;

namespace gridpp {
    static const float MV = -999;
    static const float pi = 3.14159265;
    static double radius_earth = 6.37e6;

    class KDTree;
    class Points;
    class Grid;
    std::string version();

    /** Optimal interpolation
      * @param: input 2D field of background values
      * @param: bgrid grid corresponding to input
      * @param: pobs vector of observations
      * @param: pci vector of ci values
      * @param: points observation points
    */
    int optimal_interpolation(const vec2& input,
            const gridpp::Grid& bgrid,
            const vec& pobs,
            const vec& pci,
            const gridpp::Points& points,
            float minRho,
            float hlength,
            float vlength,
            float wmin,
            int maxPoints,
            float elevGradient,
            float epsilon,
            vec2& output);

    /** Optimal interpolation for ensemble
      * @param: input 3D field of background values (Y, X, E)
      * @param: bgrid grid corresponding to input
      * @param: pobs vector of observations
      * @param: pci vector of ci values
      * @param: points observation points
    */
    int optimal_interpolation_ens(const vec3& input,
            const gridpp::Grid& bgrid,
            const vec& pobs,
            const vec& pci,
            const gridpp::Points& points,
            vec2& output);

    /** Neighbourhood in space
      * @param: input Deterministic values with dimensions Y, X
      * @param: radius filter radius in number of gridpoints
      * @param: operation one of min, mean, median, max, std
    */
    vec2 neighbourhood(const vec2& input, int radius, std::string operation);

    /** Neighbourhood in space and across members
    */
    vec2 neighbourhood_ens(const vec3& input, int radius, std::string operation);

    /** Neighbourhood quantile in space and across members
      * @param: input Ensemble values with dimensions Y, X, E
      * @param: quantile which quantile to compute (between 0 and 1)
      * @param: radius filter radius in number of gridpoints
      * @param: num_threhsolds number of thresholds to use to approximate value (0 for no approximation).
    */
    vec2 neighbourhood_quantile_ens(const vec3& input, float quantile, int radius, int num_thresholds);
    vec2 neighbourhood_quantile(const vec2& input, float quantile, int radius, int num_thresholds);

    /** Fill in values inside or outside a set of circles
      * @param: input Deterministic values with dimensions Y, X
      * @param: radiii Circle radii for each point
      * @param: value Fill in this value
      * @param: outside if True, fill outside circles, if False, fill inside circles
    */
    vec2 fill(const Grid& igrid, const vec2& input, const Points& points, const vec& radii, float value, bool outside);

    // Grid to grid interpolation
    vec2 bilinear(const Grid& igrid, const Grid& ogrid, const vec2 ivalues);
    vec2 nearest(const Grid& igrid, const Grid& ogrid, const vec2 ivalues);

    // Grid to point interpolation
    vec bilinear(const Grid& igrid, const Points& opoints, const vec2 ivalues);
    vec nearest(const Grid& igrid, const Points& opoints, const vec2 ivalues);

    namespace util {
        enum StatType {
            StatTypeMean      = 0,
            StatTypeMin       = 10,
            StatTypeMedian    = 20,
            StatTypeMax       = 30,
            StatTypeQuantile  = 40,
            StatTypeStd       = 50,
            StatTypeSum       = 60
        };
        StatType getStatType(std::string iName);
        double clock();
        void debug(std::string string);
        void error(std::string string);
        bool is_valid(float value);
        float calculate_stat(const std::vector<float>& iArray, StatType iStatType, float iQuantile=MV);
        int num_missing_values(const vec2& iArray);
        int get_lower_index(float iX, const std::vector<float>& iValues);
        int get_upper_index(float iX, const std::vector<float>& iValues);
        float interpolate(float x, const std::vector<float>& iX, const std::vector<float>& iY);

        /** Get reasonably spaced quantiles from a vector of values, ignoring duplicate values
          *  but including the first number after duplicated values. Include the lowest and highest
          *  values.
          *  @param values vector of values (unsorted, and no invalid values)
          *  @param num number of thresholds to get
        */
        vec calc_even_quantiles(const vec& values, int num);
    }

    class KDTree {
        public:
            KDTree(vec lats, vec lons);
            KDTree() {};

            /** Find single nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             * */
            int get_nearest_neighbour(float lat, float lon) const;

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            ivec get_neighbours(float lat, float lon, float radius) const;

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             * */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const;

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            int get_num_neighbours(float lat, float lon, float radius) const;

            /** Find a set of nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param num Number of points to find
             * */
            ivec get_closest_neighbours(float lat, float lon, int num) const;


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
            static float rad2deg(float deg);
            static float calc_distance(float lat1, float lon1, float lat2, float lon2);
            static float calc_distance(float x0, float y0, float z0, float x1, float y1, float z1);
            vec get_lats() const;
            vec get_lons() const;
            int size() const;
        protected:
            typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
            typedef std::pair<point, unsigned> value;
            typedef boost::geometry::model::box<point> box;
            boost::geometry::index::rtree< value, boost::geometry::index::quadratic<16> > mTree;
            vec mLats;
            vec mLons;
    };

    class Points  {
        public:
            Points(vec lats, vec lons, vec elevs=vec(), vec lafs=vec());
            int get_nearest_neighbour(float lat, float lon) const;
            ivec get_neighbours(float lat, float lon, float radius) const;
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const;
            int get_num_neighbours(float lat, float lon, float radius) const;
            ivec get_closest_neighbours(float lat, float lon, int num) const;

            vec get_lats() const;
            vec get_lons() const;
            vec get_elevs() const;
            vec get_lafs() const;
            int size() const;
        private:
            KDTree mTree;
            vec mLats;
            vec mLons;
            vec mElevs;
            vec mLafs;
    };

    class Grid {
        public:
            Grid(vec2 lats, vec2 lons, vec2 elevs=vec2(), vec2 lafs=vec2());
            ivec get_nearest_neighbour(float lat, float lon) const;
            ivec2 get_neighbours(float lat, float lon, float radius) const;
            ivec2 get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const;
            int get_num_neighbours(float lat, float lon, float radius) const;
            ivec2 get_closest_neighbours(float lat, float lon, int num) const;

            vec2 get_lats() const;
            vec2 get_lons() const;
            vec2 get_elevs() const;
            vec2 get_lafs() const;
            ivec size() const;
        private:
            KDTree mTree;
            int mX;
            vec2 get_2d(vec input) const;
            ivec get_indices(int index) const;
            ivec2 get_indices(ivec indices) const;
            vec2 mLats;
            vec2 mLons;
            vec2 mElevs;
            vec2 mLafs;
    };


    /*
    namespace downscaling {
        vec2 nearestNeighbour(const vec2& iLats, const vec2& iLons, const vec2& iValues, const vec2& oLats, const vec2& oLons, int iRadius, std::string iStatType, std::string iQuantile);
        vec2 nearestNeighbour(Field ifield, Grid ogrid, int iRadius, std::string iStatType, std::string iQuantile);
    }
    namespace calibration {
        vec2 quantileMapping(const vec2& iValues, std::string iExtrapolation, vec iQuantiles, const ParameterSet& iParameterSet);
        vec2 quantileMapping(const vec2& iValues, std::string iExtrapolation, vec iQuantiles, const Points& iParameterPoints, const vec2 iParameters);
    }
    namespace diagnose {
        vec2 computeRh(const vec2& iTemperatures, const vec2& iDewPointTemperatures);
        float computeRh(float iTemperatures, float iDewPointTemperatures);
    }

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
            int maxPoints,
            float elevGradient,
            float epsilon,
            vec2& output);

    */
};
#endif
