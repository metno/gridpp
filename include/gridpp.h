#ifndef API_H
#define API_H
#include <vector>
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
typedef std::vector<std::vector<double> > dvec2;
typedef std::vector<float> vec;
// typedef std::vector<float> fvec;
typedef std::vector<int> ivec;
typedef std::vector<std::vector<int> > ivec2;

namespace gridpp {
    // Constants
    // static const float MV = NAN;
    static const float MV = -999;
    static const float pi = 3.14159265;
    static double radius_earth = 6.37e6;

    class KDTree;
    class Points;
    class Grid;
    class Interpolator;
    class Nearest;
    class Parameters;
    /** Methods for extrapolating outside a curve */
    enum Extrapolation {
            ExtrapolationOneToOne = 0,      /**< Continue using a slope of 1 */
            ExtrapolationMeanSlope = 10,    /**< Continue outside the curve using the mean slope of the curve*/
            ExtrapolationNearestSlope = 20, /**< Continue using the slope of the last two points in the curve*/
            ExtrapolationZero = 30,         /**< Continue using a slope of 0 */
        };

    enum Operation {
        Mean      = 0,
        Min       = 10,
        Median    = 20,
        Max       = 30,
        Quantile  = 40,
        Std       = 50,
        Sum       = 60
    };
    /** The @return The gridpp version */
    std::string version();

    /****************************************
     * Data assimilation methods            *
     ****************************************/

    /** Optimal interpolation
      * @param input 2D field of background values
      * @param bgrid grid corresponding to input
      * @param pobs vector of observations
      * @param pci vector of ci values
      * @param points observation points
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
      * @param input 3D field of background values (Y, X, E)
      * @param bgrid grid corresponding to input
      * @param pobs vector of observations
      * @param pci vector of ci values
      * @param points observation points
    */
    int optimal_interpolation_ens(const vec3& input,
            const gridpp::Grid& bgrid,
            const vec& pobs,
            const vec& pci,
            const gridpp::Points& points,
            vec2& output);

    /****************************************
     * Neighbourhood methods                *
     ****************************************/

    /** Spatial neighbourhood filter
      * @param input 2D grid of values
      * @param radius Filter radius in number of gridpoints
      * @param operation One of min, mean, median, max, std
    */
    vec2 neighbourhood(const vec2& input, int radius, std::string operation);

    /** Neighbourhood filter in space and across ensemble members
      * @param input 3D vector with dimensions (Y, X, ensemble)
      * @param radius Filter radius in number of gridpoints
      * @param operation One of min, mean, median, max, std
    */
    vec2 neighbourhood_ens(const vec3& input, int radius, std::string operation);

    /** Spatial neighbourhood filter for quantile operation
      * @param input 2D grid of values
      * @param quantile Quantile to compute (between 0 and 1)
      * @param radius Filter radius in number of gridpoints
      * @param operation One of min, mean, median, max, std
    */
    vec2 neighbourhood_quantile(const vec2& input, float quantile, int radius, int num_thresholds);

    /** Neighbourhood filter space and across members for quantile operation
      * @param input 3D vector with dimensions (Y, X, ensemble)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param radius Filter radius in number of gridpoints
      * @param num_thresholds Number of thresholds to use to approximate value (0 for no approximation)
    */
    vec2 neighbourhood_quantile_ens(const vec3& input, float quantile, int radius, int num_thresholds);

    /** Neighbourhood filter space and across members for quantile operation
      * @param input 3D vector with dimensions (Y, X, ensemble)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param radius Filter radius in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_ens(const vec3& input, float quantile, int radius, const vec& thresholds);

    /** Spatial neighbourhood filter without any shortcuts. This is likely quite slow.
     *  @param input 2D grid of values
     *  @param radius Filter radius in number of gridpoints
     *  @param operation one of min, mean, median, max
    */
    vec2 neighbourhood_brute_force(const vec2& input, int radius, std::string operation);

    /** Spatial neighbourhood filter without any shortcuts for quantile operation. This is likely quite slow.
     *  @param input 2D grid of values
     *  @param radius Filter radius in number of gridpoints
     *  @param quantile Quantile to calculate for (between 0 and 1)
    */
    vec2 neighbourhood_quantile_brute_force(const vec2& input, int radius, float quantile);

    /****************************************
     * Calibration methods                  *
     ****************************************/
    // With parameter interpolator
    vec2 quantile_mapping(const Grid& grid, const vec2& input, gridpp::Extrapolation policy, const Parameters& parameters);
    // Constant map over the whole domain
    vec2 quantile_mapping(const vec2& input, const vec& x, const vec& y, gridpp::Extrapolation policy);
    // Processes a vector
    vec quantile_mapping(const vec& input, const vec& x, const vec& y, gridpp::Extrapolation policy);
    float quantile_mapping(float input, const vec& fcst, const vec& ref, gridpp::Extrapolation policy);

    /** Fill in values inside or outside a set of circles
      * @param input Deterministic values with dimensions Y, X
      * @param radiii Circle radii for each point
      * @param value Fill in this value
      * @param outside if True, fill outside circles, if False, fill inside circles
    */
    vec2 fill(const Grid& igrid, const vec2& input, const Points& points, const vec& radii, float value, bool outside);

    /****************************************
     * Downscaling methods                  *
     ****************************************/
    // Grid to grid interpolation
    vec2 nearest(const Grid& igrid, const Grid& ogrid, const vec2 ivalues);
    // Grid to point interpolation
    vec nearest(const Grid& igrid, const Points& opoints, const vec2 ivalues);

    vec2 bilinear(const Grid& igrid, const Grid& ogrid, const vec2 ivalues);
    vec bilinear(const Grid& igrid, const Points& opoints, const vec2 ivalues);
    // vec2 smart(const Grid& igrid, const Grid& ogrid, const vec2 ivalues, int radius, float num, float min_elev_diff);
    // vec2 gradient(const Grid& igrid, const Grid& ogrid, const vec2 ivalues)

    /****************************************
     * Diagnosing methods                   *
     ****************************************/

    /** Diagnose wind speed from its components
     *  @param xwind X-component of wind [any unit]
     *  @param ywind Y-component of wind [any unit]
     *  @return Wind speed [any unit]
     * */
    float wind_speed(float xwind, float ywind);
    vec wind_speed(const vec& xwind, const vec& ywind);
    vec2 wind_speed(const vec2& xwind, const vec2& ywind);

    /** Compute dewpoint temperature from temperature and relative humidity
     *  @param temperature Temperature [K]
     *  @param relative_humidity Relative humidity [1]
     *  @returns Dewpoint temperature [K]
    */
    float dewpoint(float temperature, float relative_humidity);
    vec dewpoint(const vec& temperature, const vec& relative_humidity);

    /** Compute relative humidity from temperature and dewpoint temperature
     *  @param temperature Temperature [K]
     *  @param dewpoint Dewpoint temperature [K]
     *  @returns Relative humidity [1]
    */
    float relative_humidity(float temperature, float dewpoint);

    /** Compute wetbulb temperature from temperature, pressure, and relative humidity
     *  @param temperature Temperature [K]
     *  @param pressure Air pressure [pa]
     *  @param Relative humidity [1]
     *  @returns Wetbulb temperature [K]
    */
    float wetbulb(float temperature, float pressure, float relative_humidity);

    /** Calculate pressure at a new elevation
     *  @param ielev Elevation at start point
     *  @param oelev Elevation at new point
     *  @param ipressure Pressure at start point
     *  @param itemperature Temperature at start point
     *  @return Pressure at new point
     */
    float pressure(float ielev, float oelev, float ipressure, float itemperature=288.15);

    /** Diagnose QNH from pressure and altitude
     *  @param pressure Pressure at point [pa]
     *  @param altitude Altitude of point [m]
     *  @returns QNH [pa]
    */
    float qnh(float pressure, float altitude);
    vec qnh(const vec& pressure, const vec& altitude);

    namespace util {
        // vec2 calc_gradient(const vec2& values, const vec2& aux, int radius);
        // ivec regression(const vec& x, const vec& y);
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
        void not_implemented_error();

        /** Get reasonably spaced quantiles from a vector of values, ignoring duplicate values
          *  but including the first number after duplicated values. Include the lowest and highest
          *  values.
          *  @param values vector of values (unsorted, and no invalid values)
          *  @param num number of thresholds to get
        */
        vec calc_even_quantiles(const vec& values, int num);

        /** Raises an exception if the two arrays are not the same size */
        void check_equal_size(const vec& v1, const vec& v2);
    }

    /** Helper class for Grid and points */
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

    /** Represents a vector of locations and their metadata */
    class Points  {
        public:
            Points();
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

    /** Represents a 2D grid of locations and their metadata */
    class Grid {
        public:
            Grid();
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

    /** Class for interpolating parameters at points to new points */
    class Interpolator {
        public:
            virtual float interpolate(const Points& points, const vec& values, float lat, float lon, float altitude, float land_area_fraction) const = 0;
            virtual Interpolator* clone() const = 0;
    };
    class Nearest: public Interpolator {
        public:
            Nearest(int num=1);
            float interpolate(const Points& points, const vec& values, float lat, float lon, float altitude=gridpp::MV, float land_area_fraction=gridpp::MV) const;
            Interpolator* clone() const;
        private:
            int m_num;
    };
    /** Class for encapulating parameter sets for locations and estimating them for arbitrary points */
    class Parameters {
        public:
            Parameters(const gridpp::Points& points, const vec2& values, gridpp::Interpolator& interpolator);
            ~Parameters();
            /** Estimate parameters at a given point
             *  @return Parameter set
             */
            vec get(float lat, float lon, float altitude, float land_area_fraction) const;
        private:
            Points m_points;
            vec2 m_values;
            Interpolator* m_interpolator;
    };
};
#endif
