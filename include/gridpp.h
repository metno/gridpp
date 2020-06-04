#ifndef GRIDPP_API_H
#define GRIDPP_API_H
#include <vector>
#include <string>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#ifdef _OPENMP
    #include <omp.h>
#endif

#define GRIDPP_VERSION "0.4.0a1"
#define __version__ GRIDPP_VERSION

typedef std::vector<int> ivec;      /**< 1D vector of ints */
typedef std::vector<float> vec;     /**< 1D vector of floats */
typedef std::vector<double> dvec;   /**< 1D vector of doubles */
typedef std::vector<ivec> ivec2;    /**< 2D vector of ints */
typedef std::vector<vec> vec2;      /**< 2D vector of floats */
typedef std::vector<dvec> dvec2;    /**< 2D vector of doubles */
typedef std::vector<vec2> vec3;     /**< 3D vector of floats */

namespace gridpp {
    // Constants
    static const float MV = NAN;
    static const float MV_CML = -999;
    static const float pi = 3.14159265;
    static double radius_earth = 6.37e6;

    class KDTree;
    class Points;
    class Grid;
    class Interpolator;
    class Nearest;
    class Parameters;
    class StructureFunction;
    class Transform;

    /** Methods for extrapolating outside a curve */
    enum Extrapolation {
            OneToOne = 0,      /**< Continue using a slope of 1 */
            MeanSlope = 10,    /**< Continue outside the curve using the mean slope of the curve*/
            NearestSlope = 20, /**< Continue using the slope of the last two points in the curve*/
            Zero = 30,         /**< Continue using a slope of 0 */
        };

    /** Statistical operations to reduce a vector to a scalar */
    enum Statistic {
        Mean      = 0,
        Min       = 10,
        Median    = 20,
        Max       = 30,
        Quantile  = 40,
        Std       = 50,
        Variance  = 60,
        Sum       = 70,
        Unknown   = -1
    };
    /** Convert operation string to enum */
    Statistic get_statistic(std::string name);

    /** The gridpp version
     * @return The gridpp version
    */
    std::string version();

    float* test_array(float* v, int n);

    /****************************************
     * Data assimilation methods            *
     ****************************************/

    /** Optimal interpolation for a deterministic field
      * @param bgrid Grid of background field
      * @param background 2D field of background values
      * @param points Points of observations
      * @param pobs Vector of observations
      * @param pratios Vector of ratio of observation error variance to background variance
      * @param pbackground Background with observation operator
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
    */
    vec2 optimal_interpolation(const gridpp::Grid& bgrid,
            const vec2& background,
            const gridpp::Points& points,
            const vec& pobs,
            const vec& pratios,
            const vec& pbackground,
            const gridpp::StructureFunction& structure,
            int max_points);

    vec2 optimal_interpolation_transform(const gridpp::Grid& bgrid,
            const vec2& background,
            float bsigma,
            const gridpp::Points& points,
            const vec& pobs,
            const vec& psigma,
            const vec& pbackground,
            const gridpp::StructureFunction& structure,
            int max_points,
            const gridpp::Transform& transform);

    /** Optimal interpolation using a structure function based on an ensemble 
      * See Lussana et al 2019 (DOI: 10.1002/qj.3646)
      * @param input 3D field of background values (Y, X, E)
      * @param bgrid grid corresponding to input
      * @param pobs vector of observations
      * @param pci vector of ci values
      * @param points observation points
    */
    vec3 optimal_interpolation_ensi(const gridpp::Grid& bgrid,
            const vec3& input,
            const gridpp::Points& points,
            const vec& pobs,
            const vec& psigmas,
            const vec2& pbackground,
            const gridpp::StructureFunction& structure,
            int max_points);

    /****************************************
     * Neighbourhood methods                *
     ****************************************/

    /** Spatial neighbourhood filter
      * @param input 2D grid of values
      * @param radius Filter radius in number of gridpoints
      * @param statistic Statistic to compute
    */
    vec2 neighbourhood(const vec2& input, int radius, Statistic statistic);

    /** Neighbourhood filter in space and across ensemble members
      * @param input 3D vector with dimensions (Y, X, ensemble)
      * @param radius Filter radius in number of gridpoints
      * @param statistic Statistic to compute
    */
    vec2 neighbourhood_ens(const vec3& input, int radius, Statistic statistic);

    /** Spatial neighbourhood filter. An exact but slow algorithm.
     *  @param input 2D grid of values
     *  @param radius Filter radius in number of gridpoints
     *  @param quantile Quantile to calculate for (between 0 and 1)
    */
    vec2 neighbourhood_quantile(const vec2& input, float quantile, int radius);

    /** Neighbourhood filter in space and across ensemble members. An exampt but slow algorithm.
     *  @param input 3D grid of values
     *  @param radius Filter radius in number of gridpoints
     *  @param quantile Quantile to calculate for (between 0 and 1)
    */
    vec2 neighbourhood_quantile_ens(const vec3& input, float quantile, int radius);

    /** Approximate spatial neighbourhood filter for quantile operation.
      * @param input 2D grid of values
      * @param quantile Quantile to compute (between 0 and 1)
      * @param radius Filter radius in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec2& input, float quantile, int radius, const vec& thresholds);

    /** Approximate neighbourhood filter space and across members for quantile operation
      * @param input 3D vector with dimensions (Y, X, ensemble)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param radius Filter radius in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_ens_fast(const vec3& input, float quantile, int radius, const vec& thresholds);

    /** Spatial neighbourhood filter without any shortcuts. This is quite slow and is only useful for testing.
     *  @param input 2D grid of values
     *  @param radius Filter radius in number of gridpoints
     *  @param operation one of min, mean, median, max
    */
    vec2 neighbourhood_brute_force(const vec2& input, int radius, Statistic statistic);

    /** Calculate appropriate approximation thresholds for neighbourhood quantile
     *  @param input 2D grid of values
     *  @param num_thresholds Number of thresholds
    */
    vec get_neighbourhood_thresholds(const vec2& input, int num_thresholds);

    /** Calculate appropriate approximation thresholds for neighbourhood quantile
     *  @param input 3D grid of values
     *  @param num_thresholds Number of thresholds
    */
    vec get_neighbourhood_thresholds(const vec3& input, int num_thresholds);

    /****************************************
     * Calibration methods                  *
     ****************************************/
    /** Quantile mapping of 2D grid using a constant quantile map
      * @param input Values on grid
      * @param x X-axis parameters values
      * @param y Y-axis parameters values
      * @param policy Extrapolation policy
      * @return Quantile-mapped values
    */
    vec2 quantile_mapping(const vec2& input, const vec& x, const vec& y, gridpp::Extrapolation policy);

    /** Quantile mapping of a vector of value
      * @param input Input values
      * @param x X-axis parameters values
      * @param y Y-axis parameters values
      * @param policy Extrapolation policy
      * @return Quantile-mapped output values
      *
    */
    vec quantile_mapping(const vec& input, const vec& x, const vec& y, gridpp::Extrapolation policy);

    /** Quantile mapping of a single value
      * @param input Input value
      * @param x X-axis parameters values
      * @param y Y-axis parameters values
      * @param policy Extrapolation policy
      * @return Quantile-mapped output value
      *
    */
    float quantile_mapping(float input, const vec& ref, const vec& fcst, gridpp::Extrapolation policy);

    /** Quantile mapping of 2D grid using a parameter interpolator
      * @param grid Grid
      * @param input Values on grid
      * @param policy Extrapolation policy
      * @param parameters Parameter interpolator
      * @return Quantile-mapped values
      *
    */
    vec2 quantile_mapping(const Grid& grid, const vec2& input, gridpp::Extrapolation policy, const Parameters& parameters);

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
    /** Nearest neighbour dowscaling grid to grid
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values on the output grid
    */
    vec2 nearest(const Grid& igrid, const Grid& ogrid, const vec2 ivalues);
    /** Nearest neighbour dowscaling grid to point
      * @param igrid Input grid
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values for the output points
    */
    vec nearest(const Grid& igrid, const Points& opoints, const vec2 ivalues);
    vec2 nearest(const Grid& igrid, const Points& opoints, const vec3 ivalues);

    /** Bilinear downscaling grid to grid
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values on the output grid
    */
    vec2 bilinear(const Grid& igrid, const Grid& ogrid, const vec2 ivalues);

    /** Bilinear downscaling grid to points
      * @param igrid Input grid
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values for the output points
    */
    vec bilinear(const Grid& igrid, const Points& opoints, const vec2 ivalues);

    vec2 simple_gradient(const Grid& igrid, const Grid& ogrid, const vec2 ivalues, float elev_gradient);
    vec simple_gradient(const Grid& igrid, const Points& opoints, const vec2 ivalues, float elev_gradient);

    // vec2 smart(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, int num, const StructureFunction& structure);

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

    /** Calculate dewpoint temperature from temperature and relative humidity
     *  @param temperature Temperature [K]
     *  @param relative_humidity Relative humidity [1]
     *  @returns Dewpoint temperature [K]
    */
    float dewpoint(float temperature, float relative_humidity);

    /** Vector version of dewpoint calculation
     *  @param temperature Temperatures [K]
     *  @param relative_humidity Relative humidities [1]
     *  @returns Dewpoint temperatures [K]
    */
    vec dewpoint(const vec& temperature, const vec& relative_humidity);

    /** Calculate relative humidity from temperature and dewpoint temperature
     *  @param temperature Temperature [K]
     *  @param dewpoint Dewpoint temperature [K]
     *  @returns Relative humidity [1]
    */
    float relative_humidity(float temperature, float dewpoint);

    /** Calculate wetbulb temperature from temperature, pressure, and relative humidity
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

    /** Vector version of QNH calculation
     *  @param pressure Pressures at points [pa]
     *  @param altitude Altitudes of points [m]
     *  @returns QNH [pa]
    */
    vec qnh(const vec& pressure, const vec& altitude);

    /** Set the number of OpenMP threads to use. Overwrides OMP_NUM_THREAD env variable. */
    void set_omp_threads(int num);

    /** Sets the number of OpenMP threads to 1 if OMP_NUM_THREADS undefined */
    void initialize_omp();

    /** Helper functions */
    namespace util {
        // vec2 calc_gradient(const vec2& values, const vec2& aux, int radius);
        // ivec regression(const vec& x, const vec& y);
        double clock();
        void debug(std::string string);
        void warning(std::string string);
        void error(std::string string);
        bool is_valid(float value);
        float calc_statistic(const vec& array, Statistic statistic);
        float calc_quantile(const vec& array, float quantile);
        vec calc_statistic(const vec2& array, Statistic statistic);
        vec calc_quantile(const vec2& array, float quantile=gridpp::MV);
        int num_missing_values(const vec2& iArray);

        /** Find the index in a vector that is equal or just below a value
         *  @param iX Lookup value
         *  @param iValues Lookup vector. Must be sorted.
         *  @return The index into iValues that is equal or just below iX
        */
        int get_lower_index(float iX, const std::vector<float>& iValues);

        /** Find the index in a vector that is equal or just above a value
         *  @param iX Lookup value
         *  @param iValues Lookup vector. Must be sorted.
         *  @return The index into iValues that is equal or just above iX
        */
        int get_upper_index(float iX, const std::vector<float>& iValues);

        /** Piecewise linear interpolation.o
         *  If x is outside the range of iX, then the min/max value of iY is used
         *  @param x Interpolation to this point
         *  @param iX Vector of x-axis values. Vector must be sorted.
         *  @param iY Vector of y-axis values corresponding to iX.
         *  @return Y value corresponding to x
        */
        float interpolate(float x, const std::vector<float>& iX, const std::vector<float>& iY);
        void not_implemented_error();

        /** Initialize a vector of size Y, X, with a given value */
        vec2 init_vec2(int Y, int X, float value=gridpp::MV);

        /** Get reasonably spaced quantiles from a vector of values, ignoring duplicate values
          *  but including the first number after duplicated values. Include the lowest and highest
          *  values.
          *  @param values vector of values (unsorted, and no invalid values)
          *  @param num number of thresholds to get
        */
        vec calc_even_quantiles(const vec& values, int num);

        /** Check if the grid is the same size as the 2D vector */
        bool compatible_size(const Grid& grid, const vec2& v);
        bool compatible_size(const Grid& grid, const vec3& v);
    }

    class Point {
        public:
            Point(float lat, float lon, float elev=gridpp::MV, float laf=gridpp::MV);
            float lat;
            float lon;
            float elev;
            float laf;
    };
    /** Covariance structure function */
    class StructureFunction {
        public:
            StructureFunction(float localization_distance);
            /** Correlation between two points */
            virtual float corr(const Point& p1, const Point& p2) const = 0;
            /** Maximum distance for which an observation can have an impact (localization)
              * @return Distance [m]
            */
            float localization_distance() const;
            virtual StructureFunction* clone() const = 0;
        protected:
            float barnes_rho(float dist, float length) const;
            float cressman_rho(float dist, float length) const;
            float mLocalizationDistance;
    };
    /** Simple structure function based on distance, elevation, and land area fraction */
    class BarnesStructure: public StructureFunction {
        public:
            BarnesStructure(float h, float v=0, float w=0, float hmax=gridpp::MV);
            float corr(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
        private:
            float mH;
            float mV;
            float mW;
    };

    /** Simple structure function based on distance, elevation, and land area fraction */
    class CressmanStructure: public StructureFunction {
        public:
            CressmanStructure(float h, float v=0, float w=0);
            float corr(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
        private:
            float mH;
            float mV;
            float mW;
    };
    class CrossValidation: public StructureFunction {
        public:
            CrossValidation(StructureFunction& structure, float dist=0);
            float corr(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
        private:
            StructureFunction* m_structure;
            float m_dist;
    };

    class Transform {
        public:
            virtual float forward(float value) const = 0;
            virtual float backward(float value) const = 0;
            virtual vec2 forward(const vec2& input) const;
            virtual vec2 backward(const vec2& input) const;
    };
    class Identity : public Transform {
        public:
            float forward(float value) const;
            float backward(float value) const;
    };
    class Log : public Transform {
        public:
            float forward(float value) const;
            float backward(float value) const;
    };
    class BoxCox : public Transform {
        public:
            BoxCox(float threshold);
            float forward(float value) const;
            float backward(float value) const;
        private:
            float mThreshold;
    };

    /** Helper class for Grid and Points */
    class KDTree {
        public:
            KDTree(vec lats, vec lons);
            KDTree& operator=(KDTree other);
            KDTree(const KDTree& other);
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
            static float calc_distance_fast(float lat1, float lon1, float lat2, float lon2);
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
            Points& operator=(Points other);
            Points(const Points& other);
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
            ivec get_in_domain_indices(const Grid& grid) const;
            Points get_in_domain(const Grid& grid) const;
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
