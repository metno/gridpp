#ifndef GRIDPP_API_H
#define GRIDPP_API_H
#include <vector>
#include <string>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <exception>

#define GRIDPP_VERSION "0.8.0.dev1"
#define __version__ GRIDPP_VERSION

namespace gridpp {
    /** **************************************
     * @name Short-hand notation for vectors of different dimensions sizes
     * ***************************************/ /**@{*/
    // Preferred vector types
    /** 1D float vector */
    typedef std::vector<float> vec;
    /** 2D float vector */
    typedef std::vector<vec> vec2;
    /** 3D float vector */
    typedef std::vector<vec2> vec3;
    /** 1D integer vector */
    typedef std::vector<int> ivec;
    /** 2D integer vector */
    typedef std::vector<ivec> ivec2;
    /** 3D integer vector */
    typedef std::vector<ivec2> ivec3;

    // Only use these when double is required
    /** 1D double vector */
    typedef std::vector<double> dvec;
    /** 2D double vector */
    typedef std::vector<dvec> dvec2;
    /**@}*/

    /** **************************************
     * @name Constants
     * Functions that assimilate observations onto a gridded background
     * ***************************************/ /**@{*/
    /** Missing value indicator */
    static const float MV = NAN;
    /** Missing value indicator in gridpp command-line tool */
    static const float MV_CML = -999;
    /** Mathematical constant pi */
    static const float pi = 3.14159265;
    /** Radius of the earth [m] */
    static const double radius_earth = 6.378137e6;
    /** Constant Lapse Rate moist air standard atmosphere [K/m] */
    static const float lapse_rate=0.0065;
    /** Temperature at surface in standard atmosphere [K] */
    static const float standard_surface_temperature = 288.15;
    /** Gravitational acceleration [m/s^2] */
    static const float gravit = 9.80665;
    /** Molar Mass of Dry Air [kg/mol] */
    static const float molar_mass = 0.0289644;
    /** Universal Gas Constant [kg*m^2*s^-2/(K*mol)] */
    static const float gas_constant_mol = 8.31447;
    /** Universal Gas Constant [J/(kg*K)] */
    static const float gas_constant_si = 287.05;
    /**@}*/

    class KDTree;
    class Points;
    class Grid;
    class Point;
    class Nearest;
    class StructureFunction;
    class Transform;

    /** Methods for extrapolating outside a curve */
    enum Extrapolation {
            OneToOne = 0,      /**< Continue past the end-points using a slope of 1 */
            MeanSlope = 10,    /**< Continue past the end-points using the mean slope of the curve*/
            NearestSlope = 20, /**< Continue past the end-points using the slope of the two lowermost or uppermost points in the curve */
            Zero = 30,         /**< Continue past the end-points using a slope of 0 */
            Unchanged = 40,    /**< Keep values the way they were */
        };

    /** Statistical operations to reduce a vector to a scalar */
    enum Statistic {
        Mean      = 0,  /**< Mean of values */
        Min       = 10, /**< Minimum of values */
        Median    = 20, /**< Mean of values */
        Max       = 30, /**< Maximum of values */
        Quantile  = 40, /**< A quantile from values */
        Std       = 50, /**< Standard deviation of values */
        Variance  = 60, /**< Population variance of values */
        Sum       = 70, /**< Sum of values */
        Count     = 80, /**< Count of values */
        RandomChoice = 90, /**< Randomly pick a non-nan value */
        Unknown   = -1  /**< Unknown statistic */
    };

    /** Binary verification metrics */
    enum Metric {
        Ets       = 0,  /**< Equitable threat score */
        Ts        = 1,  /**< Threat score */
        Kss       = 20, /**< Hannsen-Kuiper skill score */
        Pc        = 30, /**< Proportion correct */
        Bias      = 40, /**< Bias */
        Hss       = 50, /**< Heidke skill score */
    };

    /** Method for statistical correction */
    enum CorrectionType {
        Qq        = 0,        /**< Quantile mapping */
        Multiplicative = 10,  /**< Multiplicative */
        Additive  = 20,       /**< Additive */
    };

    /** Types of coordinates for position of points */
    enum CoordinateType {
        Geodetic = 0,      /**< Latitude and longitude */
        Cartesian = 1,     /**< X and Y */
    };

    /** Types of methods to calculate the gradient*/
    enum GradientType {
        MinMax = 0,
        LinearRegression = 10,
    };

    /** Types of simple downscaling methods */
    enum Downscaler {
        Nearest = 0,      /**< Nearest neighour downscaler */
        Bilinear = 1,     /**< Bilinear downscaler */
    };

    /** Types of comparison operators*/
    enum ComparisonOperator {
        Lt    = 0,         /**< Lower than, < */
        Leq   = 10,        /**< Lower or equal than, <= */
        Gt    = 20,        /**< Greater than, > */
        Geq   = 30,        /**< Greater or equal than, >= */
    };

    /** **************************************
     * @name Data assimilation methods
     * Functions that merge observations with a background field
     * ***************************************/ /**@{*/

    /** Optimal interpolation for a deterministic gridded field
      * @param bgrid Grid of background field
      * @param background 2D field of background values
      * @param obs_points Points of observations
      * @param obs 1D vector of observations
      * @param variance_ratios 1D vector of ratio of observation to background error variance
      * @param background_at_points 1D vector of background at observation points
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations points
      * @returns 2D vector of analised values
    */
    vec2 optimal_interpolation(const Grid& bgrid,
            const vec2& background,
            const Points& obs_points,
            const vec& obs,
            const vec& variance_ratios,
            const vec& background_at_points,
            const StructureFunction& structure,
            int max_points,
            bool allow_extrapolation=true);

    /** Optimal interpolation for a deterministic vector of points
      * @param bpoints Points of background field
      * @param background 1D field of background values
      * @param obs_points Points of observations
      * @param obs 1D vector of observations
      * @param variance_ratios 1D vector of ratio of observation to background error variance
      * @param background_at_points 1D vector of background at observation points
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations points
      * @returns 1D vector of analised values
    */
    vec optimal_interpolation(const Points& bpoints,
            const vec& background,
            const Points& obs_points,
            const vec& obs,
            const vec& variance_ratios,
            const vec& background_at_points,
            const StructureFunction& structure,
            int max_points,
            bool allow_extrapolation=true);

    /** Optimal interpolation for a deterministic gridded field including analysis variance
      * @param bgrid Grid of background field
      * @param background 2D vector of background values
      * @param bvariance 2D vector of background variances
      * @param obs_points Points of observations
      * @param obs 1D vector of observations
      * @param obs_variance 1D vector of observation variances
      * @param background_at_points 1D vector of background at observation points
      * @param bvariance_at_points 1D vector of background variance at observation points
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param analysis_variance 2D output vector of analysis variance
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
      * @returns 2D vector of analysed values
    */
    vec2 optimal_interpolation_full(const Grid& bgrid,
            const vec2& background,
            const vec2& bvariance,
            const Points& obs_points,
            const vec& obs,
            const vec& obs_variance,
            const vec& background_at_points,
            const vec& bvariance_at_points,
            const StructureFunction& structure,
            int max_points,
            vec2& analysis_variance,
            bool allow_extrapolation=true);

    /** Optimal interpolation for a deterministic vector of points including analysis variance
      * @param bpoints Points of background field
      * @param background 1D vector of background values
      * @param bvariance 1D vector of background variance
      * @param obs_points Points of observations
      * @param obs 1D vector of observations
      * @param obs_variance 1D vector of observation variances
      * @param background_at_points 1D vector of background at observation points
      * @param bvariance_at_points 1D vector of background variance at observation points
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param analysis_variance 1D output vector of analysis variance
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
      * @returns 2D vector of analised values
    */
    vec optimal_interpolation_full(const Points& bpoints,
            const vec& background,
            const vec& bvariance,
            const Points& obs_points,
            const vec& obs,
            const vec& obs_variance,
            const vec& background_at_points,
            const vec& bvariance_at_points,
            const StructureFunction& structure,
            int max_points,
            vec& analysis_variance,
            bool allow_extrapolation=true);

    /** Optimal interpolation for an ensemble gridded field
      * See Lussana et al 2019 (DOI: 10.1002/qj.3646)
      * @param bgrid Grid of background field
      * @param background 3D vector of background values (Y, X, E)
      * @param obs_points observation points
      * @param obs 1D vector of observations
      * @param obs_standard_deviations 1D vector of observation standard deviations
      * @param background_at_points 2D vector of background at observation points (L, E)
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
      * @returns 3D vector of analised values (Y, X, E)
    */
    vec3 optimal_interpolation_ensi(const Grid& bgrid,
            const vec3& background,
            const Points& obs_points,
            const vec& obs,
            const vec& obs_standard_deviations,
            const vec2& background_at_points,
            const StructureFunction& structure,
            int max_points,
            bool allow_extrapolation=true);

    /** Optimal interpolation for an ensemble point field
      * See Lussana et al 2019 (DOI: 10.1002/qj.3646)
      * @param bpoints Points of background field
      * @param background 2D vector of background values (L, E)
      * @param obs_points Observation points
      * @param obs 1D vector of observations
      * @param obs_standard_deviations 1D vector of observation standard deviations
      * @param background_at_points 2D vector of background at observation points (L, E)
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
      * @returns 2D vector of analised values (L, E)
    */
    vec2 optimal_interpolation_ensi(const Points& bpoints,
            const vec2& background,
            const Points& obs_points,
            const vec& obs,
            const vec& obs_standard_deviations,
            const vec2& background_at_points,
            const StructureFunction& structure,
            int max_points,
            bool allow_extrapolation=true);

    /** Optimal interpolation for an ensemble point field (R bindings)
      * See Lussana et al 2019 (DOI: 10.1002/qj.3646)
      * @param bpoints Points of background field
      * @param background 2D vector of background values (L, E)
      * @param obs_points Observation points
      * @param obs 1D vector of observations
      * @param obs_standard_deviations 1D vector of observation standard deviations
      * @param background_at_points 2D vector of background at observation points (L, E)
      * @param which_structfun structure function to use (0=Barnes;1=MixA)
      * @param dh length scale for the horizontal structure function
      * @param dz length scale for the vertical structure function
      * @param dw minimum value of the correlation coefficient for laf structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
      * @returns 2D vector of analised values (L, E)
    */
    vec2 R_optimal_interpolation_ensi(const Points& bpoints,
            const vec2& background,
            const Points& obs_points,
            const vec& obs,
            const vec& obs_standard_deviations,
            const vec2& background_at_points,
/*            const StructureFunction& structure, */
            int which_structfun,
            float dh,
            float dz,
            float dw,
            int max_points,
            bool allow_extrapolation=true);

    /** Optimal interpolation for an ensemble gridded field (alternative version)
      * Work in progress
      * @param bgrid Grid of background field
      * @param background 3D vector of (left) background values to update (Y, X, E)
      * @param background 3D vector of (LEFT) background values (Y, X, E) used to compute correlations
      * @param obs_points observation points
      * @param obs 2D vector of perturbed observations (S, E)
      * @param background 3D vector of (right) background values used to compute innovations (Y, X, E)
      * @param background 3D vector of (RIGHT) background values (Y, X, E) used to compute correlations
      * @param structure Structure function
      * @param variance_ratio (ratio of observation to right background error variance)
      * @param standard deviation ratio (ratio of left to right background error standard deviation)
      * @param weight given to the analysis increment
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
      * @returns 3D vector of analised values (Y, X, E)
    */
/*    vec3 optimal_interpolation_ensi_lr(const Grid& bgrid,
            const vec3& background_l,
            const vec3& background_L,
            const Points& obs_points,
            const vec2& obs,
            const vec2& pbackground_r,
            const vec2& pbackground_R,
            const StructureFunction& structure,
            float var_ratios_or,
            float std_ratios_lr,
            float weigth,
            int max_points,
            bool allow_extrapolation=true); */

    /** Optimal interpolation for an ensemble gridded field (alternative version that works with R bindings)
      * Work in progress
      * @param bpoints Points of background field
      * @param background 2D vector of (left) background values to update (M, E) M=num. grid points
      * @param background 2D vector of (LEFT) background values (M, E) used to compute correlations
      * @param obs_points Observation points
      * @param obs 2D vector of perturbed observations (S, E) S=num. obs points
      * @param background 2D vector of (right) background values used to compute innovations (S, E)
      * @param background 2D vector of (RIGHT) background values (S, E) used to compute correlations
      * @param which_structfun structure function to use (0=Barnes;1=MixA)
      * @param dh length scale for the horizontal structure function
      * @param dz length scale for the vertical structure function
      * @param dw minimum value of the correlation coefficient for laf structure function
      * @param variance_ratio (ratio of observation to right background error variance)
      * @param standard deviation ratio (ratio of left to right background error standard deviation)
      * @param weight given to the analysis increment
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
      * @returns 2D vector of analised values (M, E)
    */
    vec2 R_optimal_interpolation_ensi_lr(const Points& bpoints,
            const vec2& background_l,
            const vec2& background_L,
            const Points& obs_points,
            const vec2& obs,
            const vec2& pbackground_r,
            const vec2& pbackground_R,
/*            const StructureFunction& structure, */
            int which_structfun,
            float dh,
            float dz,
            float dw,
            float var_ratios_or,
            float std_ratios_lr,
            float weigth,
            int max_points,
            bool allow_extrapolation=true);

    /** Correction of a gridded field ensuring the distribution of values nearby match that of
      * observations. This is an experimental method.
      * @param bgrid Grid of background field
      * @param background 2D vector of background values (Y, X)
      * @param obs_points Points of observations
      * @param obs 1D vector of observations
      * @param background_at_points 1D vector of background values at observation points
      * @param structure Structure function
      * @param min_quantile Truncate quantile map below this quantile level
      * @param max_quantile Truncate quantile map above this quantile level
      * @param max_points Maximum number of points used within localization radius (not used at moment)
      * @returns 2D vector of corrected values
    */
    vec2 local_distribution_correction(const Grid& bgrid,
            const vec2& background,
            const Points& obs_points,
            const vec& obs,
            const vec& background_at_points,
            const StructureFunction& structure,
            float min_quantile,
            float max_quantile,
            int min_points=0);

    /** Correction of a gridded field ensuring the distribution of values nearby match that of
      * observations. This is an experimental method. Version with multiple number of timesteps. Pool
      * all observations across time in when computing the calibration curve.
      * @param bgrid Grid of background
      * @param background 2D vector of background values (Y, X)
      * @param obs_points Points of observations
      * @param obs 2D vector of observations (Time, Points)
      * @param background_at_points 2D vector of background values at observation points (Time, Points)
      * @param structure structure function
      * @param min_quantile Truncate quantile map below this quantile
      * @param max_quantile Truncate quantile map above this quantile
      * @param max_points Maximum number of points used within localization radius (not used at moment)
      * @returns 2D vector of corrected values
    */
    vec2 local_distribution_correction(const Grid& bgrid,
            const vec2& background,
            const Points& obs_points,
            const vec2& obs,
            const vec2& background_at_points,
            const StructureFunction& structure,
            float min_quantile,
            float max_quantile,
            int min_points=0);

    /** Fill in values inside or outside a set of circles (useful for masking)
      * @param igrid Grid of input
      * @param input 2D vector of values (Y, X)
      * @param points Points to fill in
      * @param radii 1D vector of circle radii for each point [m]
      * @param value Value to fill in
      * @param outside If True, fill outside circles, if False, fill inside circles
      * @returns 2D vector of final field
    */
    vec2 fill(const Grid& igrid,
            const vec2& input,
            const Points& points,
            const vec& radii,
            float value,
            bool outside);

    /** Insert observations into gridded field using a square box
      * @param grid Grid of background
      * @param background 2D vector of background values (Y, X)
      * @param points Points of observations
      * @param observations 1D vector of observations
      * @param halfwidths Half width of square (in number of grid points) where observations are inserted for each point
      * @param max_elev_diff Only insert where elevation difference between grid and point is less than this value [m]
      * @returns 2D vector of final field
    */
    vec2 doping_square(const Grid& igrid,
            const vec2& background,
            const Points& points,
            const vec& observations,
            const ivec& halfwidths,
            float max_elev_diff=gridpp::MV);

    /** Insert observations into gridded field using a circle
      * @param grid Grid of background
      * @param background 2D vector of background values with (Y, X)
      * @param points Points of observations
      * @param observations 1D vector of observations
      * @param radii Radius of circle where observations are inserted for each point [m]
      * @param max_elev_diff Only insert where elevation difference between grid and point is less than this value
      * @returns 2D vector of final field
    */
    vec2 doping_circle(const Grid& igrid,
            const vec2& background,
            const Points& points,
            const vec& observations,
            const vec& radii,
            float max_elev_diff=gridpp::MV);

    /** **************************************
     * @name Distributions
     * Functions that extract values from probability distributions
     * ***************************************/ /**@{*/

    /** Extract quantiles from a gamma distribution
     * @param levels 1D vector of quantile levels to retrieve
     * @param shape 1D vector of shape parameter of gamma distribution
     * @param scale 1D vector of scale parameter of gamma distribution
     * @returns 1D vector of quantiles
    */
    vec gamma_inv(const vec& levels, const vec& shape, const vec& scale);

    /**@}*/

    /** **************************************
     * @name Spatial neighbourhood filters
     * Functions that apply neighbourhood filters on a gridded field
     * ***************************************/ /**@{*/

    /** Spatial neighbourhood filter, computing a statistic for a sliding square window
      * @param input 2D vector of input values (Y, X)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param statistic Statistic to compute
      * @returns 2D vector of final field (Y, x)
    */
    vec2 neighbourhood(const vec2& input, int halfwidth, Statistic statistic);

    /** Spatial neighbourhood filter for an ensemble of fields. Computes static across neighbourhood
      * and all members
      * @param input 3D vector of input values (Y, X, E)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param statistic Statistic to compute
      * @returns 2D vector of final field (Y, x)
    */
    vec2 neighbourhood(const vec3& input, int halfwidth, Statistic statistic);

    /** Computes a quantile in a sliding square neighbourhood
      * @param input 2D vector of values
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @returns 2D vector of final field (Y, x)
    */
    vec2 neighbourhood_quantile(const vec2& input, float quantile, int halfwidth);

    /** Computes a quantile in a sliding square neighbourhood across ensemble members
      * @param input 3D vector of values (Y, X, E)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @returns 2D vector of final field (Y, x)
    */
    vec2 neighbourhood_quantile(const vec3& input, float quantile, int halfwidth);

    /** Fast and approximate neighbourhood quantile
      * @param input 2D vector of values (Y, X)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds 1D vector of thresholds to use to approximate value
      * @returns 2D vector of final field (Y, X)
    */
    vec2 neighbourhood_quantile_fast(const vec2& input, float quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile for ensemble of fields
      * @param input 3D vector of values (Y, X, E)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds 1D vector of thresholds to use to approximate value
      * @returns 2D vector of final field (Y, X)
    */
    vec2 neighbourhood_quantile_fast(const vec3& input, float quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile, with spatially varying quantile levels
      * @param input 2D vector of values (Y, X)
      * @param quantile 2D vector quantile levels to compute (between 0 and 1) (Y, X)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds 1D vector of thresholds to use to approximate value
      * @returns 2D vector of final field (Y, X)
    */
    vec2 neighbourhood_quantile_fast(const vec2& input, const vec2& quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile for ensemble of fields, with spatially varying quantile
      * @param input 3D vector of values with dimensions (Y, X, E)
      * @param quantile 2D vector of quantile levels to compute (between 0 and 1) (Y, X)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds 1D vector of thresholds to use to approximate value
      * @returns 2D vector of final field (Y, X)
    */
    vec2 neighbourhood_quantile_fast(const vec3& input, const vec2& quantile, int halfwidth, const vec& thresholds);

    /** Spatial neighbourhood filter without any shortcuts. This is quite slow and is only useful for testing.
      * @param input 2D vector of values (Y, X)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param statistic Statistic to compute
      * @returns 2D vector of final field (Y, X)
    */
    vec2 neighbourhood_brute_force(const vec2& input, int halfwidth, Statistic statistic);

    /** Spatial ensemble neighbourhood filter without any shortcuts. This is quite slow and is only useful for testing.
      * @param input 3D vector of values (Y, X, E)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param statistic Statistic to compute
      * @returns 2D vector of final field (Y, X)
    */
    vec2 neighbourhood_brute_force(const vec3& input, int halfwidth, Statistic statistic);

    /** Calculate appropriate approximation thresholds for fast neighbourhood quantile
      * @param input 2D vector of values (Y, X)
      * @param num_thresholds Number of thresholds
      * @returns 1D vector of approximation thresholds
    */
    vec get_neighbourhood_thresholds(const vec2& input, int num_thresholds);

    /** Calculate appropriate approximation thresholds for neighbourhood quantile based on an ensemble
      * @param input 3D vector of values (Y, X, T)
      * @param num_thresholds Number of thresholds
      * @returns 1D vector of approximation thresholds
    */
    vec get_neighbourhood_thresholds(const vec3& input, int num_thresholds);

    /** Computes gradients based on values in neighbourhood
     *  @param base 2D vector of dependent values, missing values are ignored (Y, X)
     *  @param values 2D vector of independent values, missing values are ignored (Y, X)
     *  @param gradient_type What gradient type to compute
     *  @param halfwidth Neighbourhood half width in number of gridpoints
     *  @param min_nim Minimum number of points required to compute gradient
     *  @param min_range Minimum range of base to compute gradient
     *  @param default_gradient Use this gradient if minimum number of points not available
     *  @returns 2D vector of gradients (Y, X)
    */
    vec2 calc_gradient(const vec2& base, const vec2& values, GradientType gradient_type, int halfwidth, int min_num=2, float min_range=gridpp::MV, float default_gradient=0);

    /** Find suitable value in neighbourhood that match a search criteria. If values are found
     * within the range, then the average of the values are used. If not, the closest point to the
     * range is used, provided this is within search_delta of the criteria range. if no point
     * fulfills this, the original value is not modified.
     * @param input 2D vector of input values, e.g. temperatures (Y, X)
     * @param search_array 2D vector to search in, e.g. land area fraction (Y, X)
     * @param halfwidth Neighbourhood halfwidth to search in
     * @param search_target_min Lower bound of criteria, e.g. a land area fraction value
     * @param search_target_Max Upper bound of criteria, e.g. a land area fraction value
     * @param search_delta Max distance outside target range where the most suitable point is still used
     * @param apply_array 2D integer vector of 1 (correct this grid point) or 0 (don't correct)
     * @returns 2D vector of output values (Y, X)
     */
    vec2 neighbourhood_search(const vec2& array, const vec2& search_array, int halfwidth, float search_target_min, float search_target_max, float search_delta, const ivec2& apply_array=ivec2());

    /** Deprecated: Compute neighbourhood statistic on ensemble field
      * @deprecated Use neighbourhood() function */
    vec2 neighbourhood_ens(const vec3& input, int halfwidth, Statistic statistic);
    /** Deprecated: Compute neighbourhood quantiles on ensemble field
      * @deprecated Use neighbourhood_quantile() function */
    vec2 neighbourhood_quantile_ens(const vec3& input, float quantile, int halfwidth);
    /** Deprecated: Compute neighbourhood quantiles fast on ensemble field
      * @deprecated Use neighbourhood_quantile_fast() function */
    vec2 neighbourhood_quantile_ens_fast(const vec3& input, float quantile, int radius, const vec& thresholds);
    /**@}*/

    /** **************************************
     * @name Calibration methods              
     * Functions that apply statistical methods to data
     * ***************************************/ /**@{*/

    /** Create quantile mapping calibration curve
     *  @param ref 1D vector of reference values (observations)
     *  @param fcst 1D vector of forecast values
     *  @param output_fcst 1D vector of output forecast quantiles
     *  @param quantiles 1D vector of quantile levels to extract. If empty, use all values.
     *  @returns 1D vector of output reference quantiles
    */
    vec quantile_mapping_curve(const vec& ref, const vec& fcst, vec& output_fcst, vec quantiles=vec());

    /** Create calibration curve that optimizes a metric
     *  @param ref 1D vector of reference values (observations)
     *  @param fcst 1d vector of forecast values
     *  @param thresholds 1D vector of thresholds for computing optimal values for
     *  @param metric A Metric to optimize for
     *  @param output_fcst 1D vector of output forecast quantiles
     *  @returns 1D vector of output reference quantiles
    */
    vec metric_optimizer_curve(const vec& ref, const vec& fcst, const vec& thresholds, Metric metric, vec& output_fcst);

    /** Apply arbitrary calibration curve to a single value
     *  @param fcst 1D vector of forecast values to correct
     *  @param curve_ref 1D vector of reference curve values
     *  @param curve_fcst 1D vector of forecast curve values
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @returns Calibrated forecasts
    */
    float apply_curve(float fcst, const vec& curve_ref, const vec& curve_fcst, Extrapolation policy_below, Extrapolation policy_above);

    /** Apply arbitrary calibration curve to 1D forecasts
     *  @param fcst 1D vector of forecast values to correct
     *  @param curve_ref 1D vector of reference curve values
     *  @param curve_fcst 1D vector of forecast curve values
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @returns 1D vector of calibrated forecasts
    */
    vec apply_curve(const vec& fcst, const vec& curve_ref, const vec& curve_fcst, Extrapolation policy_below, Extrapolation policy_above);

    /** Apply arbitrary calibration curve to 2D forecasts
     *  @param fcst 2D vector of forecast values to correct
     *  @param curve_ref 1D vector of reference curve values
     *  @param curve_fcst 1D vector of forecast curve values
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @returns 2D array of calibrated forecasts (Y, X)
    */
    vec2 apply_curve(const vec2& fcst, const vec& curve_ref, const vec& curve_fcst, Extrapolation policy_below, Extrapolation policy_above);

    /** Apply arbitrary calibration curve to 2D forecasts with spatially varying QQ map
     *  @param fcst 2D vector of forecast values to correct
     *  @param curve_ref 3D vector of reference curve values (Y, X, Curve)
     *  @param curve_fcst 3D vector of Forecast curve values (Y, X, Curve)
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @returns 2D vector of calibrated forecasts (Y, X)
    */
    vec2 apply_curve(const vec2& fcst, const vec3& curve_ref, const vec3& curve_fcst, Extrapolation policy_below, Extrapolation policy_above);

    /** Ensure calibration curve is monotonic, by removing points on the curve
     *  @param curve_ref 1D vector of reference curve values
     *  @param curve_fcst 1D vector of forecast curve values
     *  @param output_fcst 1D vector of output forecast curve values
     *  @returns 1D vector of output reference curve values
    */
    vec monotonize_curve(vec curve_ref, vec curve_fcst, vec& output_fcst);

    /** Compute the optimal threshold to remap an input threshold to
     *  @param curve_ref 1D vector of reference curve values
     *  @param curve_fcst 1D vector of reference curve values
     *  @param threshold Input threshold
     *  @param metric Metric to optimize
     *  @returns Optimal remapping threshold
    */
    float get_optimal_threshold(const vec& curve_ref, const vec& curve_fcst, float threshold, Metric metric);

    /** Compute the score for a 2x2 contingency table
     *  @param a Fraction of hits
     *  @param b Fraction of false alarms
     *  @param c Fraction of misses
     *  @param d Fraction of correct rejections
     *  @param metric Metric to conpute score for
     *  @returns Score
    */
    float calc_score(float a, float b, float c, float d, Metric metric);

    /** Compute the score for a 2x2 contingency table by thresholding
     *  @param ref 1D vector of reference values
     *  @param fcst 1D vector forecast values
     *  @param threshold Threshold to compute score for
     *  @param metric Metric to conpute score for
     *  @returns Score
    */
    float calc_score(const vec& ref, const vec& fcst, float threshold, Metric metric);

    /** Compute the score for a 2x2 contingency table by thresholding and allowing a different
     *  threshold for the forecast 
     *  @param ref 1D vector of reference values
     *  @param fcst 1D vector forecast values
     *  @param threshold Threshold for the reference values
     *  @param fthreshold Threshold for the forecast values
     *  @param metric Metric to conpute score for
     *  @returns Score
    */
    float calc_score(const vec& ref, const vec& fcst, float threshold, float fthreshold, Metric metric);

    /**@}*/

    /** **************************************
     * @name Downscaling methods
     * Functions that interpolate data from one grid to another
     * ***************************************/ /**@{*/

    /** Downscale a gridded field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 2D vector of values on the input grid (Yi, Xi)
      * @param downscaler Downscaling method
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 downscaling(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, Downscaler downscaler);

    /** Downscale a gridded ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 3D vector of values on the input grid (E, Yi, Xi)
      * @param downscaler Downscaling method
      * @returns 3D vector of output values (E, Yo, Xo)
    */
    vec3 downscaling(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, Downscaler downscaler);

    /** Downscale a gridded field to points
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 2D vector of values on the input grid (Yi, Xi)
      * @param downscaler Downscaling method
      * @returns 1D vector of output values
    */
    vec downscaling(const Grid& igrid, const Points& opoints, const vec2& ivalues, Downscaler downscaler);

    /** Downscale a gridded ensemble field to points
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 3D vector of values on the input grid (E, Yi, Xi)
      * @param downscaler Downscaling method
      * @returns 2D vector of output values (E, Points)
    */
    vec2 downscaling(const Grid& igrid, const Points& opoints, const vec3& ivalues, Downscaler downscaler);

    /** Nearest neighbour dowscaling grid to grid for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 2D vector of input values (Yi, Xi)
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 nearest(const Grid& igrid, const Grid& ogrid, const vec2& ivalues);

    /** Nearest neighbour dowscaling grid to grid for an ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 3D vector of input values (E, Yi, Xi)
      * @returns 3D vector of output values (E, Yo, Xo)
    */
    vec3 nearest(const Grid& igrid, const Grid& ogrid, const vec3& ivalues);

    /** Nearest neighbour dowscaling grid to point for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 2D vector of input values grid (Yi, Xi)
      * @returns 1D vector of output values
    */
    vec nearest(const Grid& igrid, const Points& opoints, const vec2& ivalues);

    /** Nearest neighbour dowscaling grid to point for an ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 3D vector of input values (E, Yi, Xi)
      * @returns 2D vector of output values (E, Points)
    */
    vec2 nearest(const Grid& igrid, const Points& opoints, const vec3& ivalues);

    /** Nearest neighbour dowscaling point to point for a deterministic field
      * @param ipoints Points of input
      * @param opoints Points of downscaled output
      * @param ivalues 1D vector of input values
      * @returns Values 1D vector of output values
    */
    vec nearest(const Points& ipoints, const Points& opoints, const vec& ivalues);

    /** Nearest neighbour dowscaling point to point for an ensemble field
      * @param ipoints Points of input (Pi)
      * @param opoints Points of downscaled output Po)
      * @param ivalues 2D vector of input values (E, Pi)
      * @returns 2D vector of output values (E, Po)
    */
    vec2 nearest(const Points& ipoints, const Points& opoints, const vec2& ivalues);

    /** Nearest neighbour dowscaling point to grid for a deterministic field
      * @param ipoints Points of input
      * @param ogrid Grid of downscaled output
      * @param ivalues 1D vector of values for the input points
      * @returns 2D vector of output values (Y, X)
    */
    vec2 nearest(const Points& ipoints, const Grid& ogrid, const vec& ivalues);

    /** Nearest neighbour dowscaling point to grid for an ensemble field
      * @param ipoints Points of input
      * @param ogrid Grid of downscaled output
      * @param ivalues 2D vector input values (E, Points)
      * @returns 3D vector of output values (E, Y, X)
    */
    vec3 nearest(const Points& ipoints, const Grid& ogrid, const vec2& ivalues);

    /** Nearest neighbour downscaling grid to grid and probability in one
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 3D vector of input values (Yi, Xi, E)
      * @param threshold 2D vector of threshold values (Yo, Xo)
      * @param comparison_operator Comparison operator to create probability
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 downscale_probability(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec2& threshold, const ComparisonOperator& comparison_operator);

    /** Masked ensemble downscaling and consensus. This method computes a high-resolution consensus
      * from a low-resolution ensemble, without downscaling all members. See the wiki for more
      * details.
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues_true 3D vector of input values to use when condition is true (Yi, Xi, E)
      * @param ivalues_false 3D vector of input values to use when condition is false (Yi, Xi, E)
      * @param threshold_values 3D vector of values (Yi, Xi, E) defining the mask for ivalues after applying the theshold
      * @param threshold 2D vector of thresholds on output grid that is upscaled to input grid (Yo, Xo)
      * @param comparison_operator Comparison operator to create probability
      * @param statistic statistic to compute over the ensemble dimension after downscaling
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 mask_threshold_downscale_consensus(const Grid& igrid, const Grid& ogrid, const vec3& ivalues_true, const vec3& ivalues_false, const vec3& theshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const Statistic& statistic);

    /** Masked ensemble downscaling and quantile extraction. This method extractes a quantile for a
      * high-resolution grid from a low-resolution ensemble, without downscaling all members. See
      * the wiki for more details.
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues_true 3D vector of input values to use when condition is true (Yi, Xi, E)
      * @param ivalues_false 3D vector of input values to use when condition is false (Yi, Xi, E)
      * @param threshold_values 3D vector of values (Yi, Xi, E) defining the mask for ivalues after applying the theshold
      * @param threshold 2D vector of thresholds on output grid that is upscaled to input grid (Yo, Xo)
      * @param comparison_operator Comparison operator to create probability
      * @param quantile_level Quantile level (between 0 and 1) to extract
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 mask_threshold_downscale_quantile(const Grid& igrid, const Grid& ogrid, const vec3& ivalues_true, const vec3& ivalues_false, const vec3& theshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const float quantile_level);

    /** Bilinear downscaling grid to grid for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 2D vector of input values (Yi, Xi)
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 bilinear(const Grid& igrid, const Grid& ogrid, const vec2& ivalues);

    /** Bilinear dowscaling grid to grid for an ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 3D vector of input values (E, Yi, Xi)
      * @returns 3D vector of output values (E, Yo, Xo)
    */
    vec3 bilinear(const Grid& igrid, const Grid& ogrid, const vec3& ivalues);

    /** Bilinear dowscaling grid to point for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 2D vector of input values grid (Yi, Xi)
      * @returns 1D vector of output values
    */
    vec bilinear(const Grid& igrid, const Points& opoints, const vec2& ivalues);

    /** Bilinear dowscaling grid to point for an ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 3D vector of input values (E, Yi, Xi)
      * @returns 2D vector of output values (E, Points)
    */
    vec2 bilinear(const Grid& igrid, const Points& opoints, const vec3& ivalues);

    /** Elevation correction dowscaling grid to grid for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 2D vector of input values (Yi, Xi)
      * @param elev_gradient Elevation gradient [input unit/m]
      * @param downscaler Downscaling method
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 simple_gradient(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, float elev_gradient, Downscaler downscaler=Nearest);

    /** Elevation correction dowscaling grid to grid for an ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 3D vector of input values (E, Yi, Xi)
      * @param elev_gradient Elevation gradient [input unit/m]
      * @param downscaler Downscaling method
      * @returns 3D vector of output values (E, Yo, Xo)
    */
    vec3 simple_gradient(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, float elev_gradient, Downscaler downscaler=Nearest);

    /** Elevation correction dowscaling grid to point for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 2D vector of input values grid (Yi, Xi)
      * @param elev_gradient Elevation gradient [input unit/m]
      * @param downscaler Downscaling method
      * @returns 1D vector of output values
    */
    vec simple_gradient(const Grid& igrid, const Points& opoints, const vec2& ivalues, float elev_gradient, Downscaler downscaler=Nearest);

    /** Elevation correction dowscaling grid to point for an ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 3D vector of input values (E, Yi, Xi)
      * @param elev_gradient Elevation gradient [input unit/m]
      * @param downscaler Downscaling method
      * @returns 2D vector of output values (E, Points)
    */
    vec2 simple_gradient(const Grid& igrid, const Points& opoints, const vec3& ivalues, float elev_gradient, Downscaler downscaler=Nearest);

    /** Compute Downscale
    *@param igrid input grid
    *@param ogrid output grid
    *@param ivalues values from igrid
    *@param elev_gradient elevation gradient
    *@param laf_gradient land area fraction gradient
    */
    /** Spatially varying gradient dowscaling grid to grid for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 2D vector of input values (Yi, Xi)
      * @param elev_gradient 2D vector of elevation gradients (Yi, Xi) [input unit/m]
      * @param laf_gradient 2D vector of land-area-fraction gradients (Yi, Xi) [input unit/1]
      * @param downscaler Downscaling method
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 full_gradient(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient=vec2(), Downscaler downscaler=Nearest);

    /** Spatially varying gradient dowscaling grid to grid for an ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 3D vector of input values (E, Yi, Xi)
      * @param elev_gradient 3D vector of elevation gradients (Yi, Xi, E) [input unit/m]
      * @param laf_gradient 3D vector of land-area-fraction gradients (Yi, Xi, E) [input unit/1]
      * @param downscaler Downscaling method
      * @returns 3D vector of output values (E, Yo, Xo)
    */
    vec3 full_gradient(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec3& elev_gradient, const vec3& laf_gradient, Downscaler downscaler=Nearest);

    /** Spatially varying gradient dowscaling grid to point for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 2D vector of input values grid (Yi, Xi)
      * @param elev_gradient 2D vector of elevation gradients (Yi, Xi) [input unit/m]
      * @param laf_gradient 2D vector of land-area-fraction gradients (Yi, Xi) [input unit/1]
      * @param downscaler Downscaling method
      * @returns 1D vector of output values
    */
    vec full_gradient(const Grid& igrid, const Points& opoints, const vec2& ivalues, const vec2& elev_gradient, const vec2& laf_gradient, Downscaler downscaler=Nearest);

    /** Spatially varying gradient dowscaling grid to point for an ensemble field
      * @param igrid Grid of input (Yi, Xi)
      * @param opoints Points of downscaled output
      * @param ivalues 3D vector of input values (E, Yi, Xi)
      * @param elev_gradient 3D vector of elevation gradients (E, Yi, Xi) [input unit/m]
      * @param laf_gradient 3D vector of land-area-fraction gradients (E, Yi, Xi) [input unit/1]
      * @param downscaler Downscaling method
      * @returns 2D vector of output values (E, Points)
    */
    vec2 full_gradient(const Grid& igrid, const Points& opoints, const vec3& ivalues, const vec3& elev_gradient, const vec3& laf_gradient, Downscaler downscaler=Nearest);

    /* Elevation and land area fraction downscaling with debug output fields. For debugging only,
    */
    vec3 full_gradient_debug(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient=vec2(), Downscaler downscaler=Nearest);

    /** Smart neighbour downscaling grid to grid for a deterministic field
      * @param igrid Grid of input (Yi, Xi)
      * @param ogrid Grid of downscaled output (Yo, Xo)
      * @param ivalues 2D vector of input values (Yi, Xi)
      * @param num Number of neighbours to average
      * @param structure Structure function for determining how similar neighbour gridpoints are
      * @returns 2D vector of output values (Yo, Xo)
    */
    vec2 smart(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, int num, const StructureFunction& structure);
    /**@}*/

    /** **************************************
     * @name Grid calculations
     * Functions that calculate statistics on a grid
     * ***************************************/ /**@{*/

    /** For each point, counts the number of gridpoints within the radius
     *  @param grid Input grid
     *  @param points Output points
     *  @param radius Radius [m]
     *  @returns 1D vector of number of gridpoints
    */
    vec count(const Grid& grid, const Points& points, float radius);

    /** For each gridpoint, counts the number of gridpoints within the radius
     *  @param igrid Input grid
     *  @param ogrid Output grid
     *  @param radius Radius [m]
     *  @returns 2D vector of number of gridpoints
    */
    vec2 count(const Grid& igrid, const Grid& ogrid, float radius);

    /** For each gridpoint, counts the number of points within the radius
     *  @param points Input points
     *  @param grid Output grid
     *  @param radius Radius [m]
     *  @returns 2D vector of number of points
    */
    vec2 count(const Points& points, const Grid& grid, float radius);

    /** For each point, counts the number of points within the radius
     *  @param ipoints Input points (Pi)
     *  @param opoints Output points (Po)
     *  @param radius Radius [m]
     *  @returns 1D vector of number of points (Po)
    */
    vec count(const Points& ipoints, const Points& opoints, float radius);

    /** For each point, calculates the distance to nearest gridpoint
     *  @param grid Grid
     *  @param points Points
     *  @param num Number of points
     *  @returns 1D vector of distance [m] to nearest gridpoint for each point
    */
    vec distance(const Grid& grid, const Points& points, int num=1);

    /** For each output gridpoint, calculate the distance to nearest input gridpoint
     *  @param grid Grid (Yi, Xi)
     *  @param ogrid Output grid /Yo, Xo)
     *  @param num Number of points
     *  @returns 2D vector of distance [m] to nearest gridpoint for each gridpoint (Yo, Xo)
    */
    vec2 distance(const Grid& igrid, const Grid& ogrid, int num=1);

    /** For each gridpoint, calculates the distance to nearest point
     *  @param points Points
     *  @param grid Grid
     *  @param num Number of points
     *  @returns 2D vector of distance [m] to nearest point for each gridpoint
    */
    vec2 distance(const Points& points, const Grid& grid, int num=1);

    /** For each output point, calculates the distance to nearest input point
     *  @param ipoints Input points (Pi)
     *  @param opoints Output points (Po)
     *  @param num Number of points
     *  @returns 1D vector of distance [m] to nearest point for each point (Po)
    */
    vec distance(const Points& ipoints, const Points& opoint, int num=1);

    /** Fill in missing values based on nearby values
      * @param values 2D vector of values
      * @returns 2D vector of values without any missing values
    */
    vec2 fill_missing(const vec2& values);

    /** Aggregate points onto a grid. Writes gridpp::MV where there are not enough points.
      * @param grid Output grid to aggregate to
      * @param points Input points
      * @param values 1D vector of input values
      * @param radius Circle radius for aggregating points [m]
      * @param min_num Minimum number of points in radius to create a value
      * @param statistic Statistic to compute on points within radius
      * @returns 2D vector of output values
    */
    vec2 gridding(const Grid& grid, const Points& points, const vec& values, float radius, int min_num, Statistic statistic);

    /** Aggregate points onto a points. Writes gridpp::MV where there are not enough points.
      * @param opoints Output points to aggregate to
      * @param ipoints Input points
      * @param values 1D vector of input Values
      * @param radius Circle radius for aggregate points [m]
      * @param min_num Minimum number of points in radius to create a value
      * @param statistic Statistic to compute on points within radius
      * @returns 1D vector of output values
    */
    vec gridding(const Points& opoints, const Points& ipoints, const vec& values, float radius, int min_num, Statistic statistic);

    /** Assign each point to nearest neighbour in grid and aggregate values. Writes gridpp::MV where
      * there are not enough observations.
      * @param grid Output grid to aggregate to
      * @param points Input points
      * @param values Input values
      * @param min_num Minimum number of points in a gridpoint to create a value
      * @param statistic Statistic to compute on points within gridpoint
      * @returns 2D vector of output values
    */
    vec2 gridding_nearest(const Grid& grid, const Points& points, const vec& values, int min_num, gridpp::Statistic statistic);

    /** Assign each ipoint to nearest neighbour in opoints and aggregate values. Writes gridpp::MV
      * where there are not enough observations.
      * @param opoints Output points to aggregate to
      * @param ipoints Input points
      * @param values Input values
      * @param min_num Minimum number of points in a gridpoint to create a value
      * @param statistic Statistic to compute on points within gridpoint
      * @returns 1D vector of output values
    */
    vec gridding_nearest(const Points& opoints, const Points& ipoints, const vec& values, int min_num, gridpp::Statistic statistic);

    /** Compute a score for a metric of all points within a radius. This is an experimental method. */
    vec2 neighbourhood_score(const Grid& grid, const Points& points, const vec2& fcst, const vec& ref, int half_width, gridpp::Metric metric, float threshold);

    /**@}*/

    /** ****************************************
     * @name Diagnosing meteorological variables
     * Functions that diagnose a meteorological variable based on other variables
     * *****************************************/ /**@{*/

    /** Calculate dewpoint temperature from temperature and relative humidity
     *  @param temperature Temperature [K]
     *  @param relative_humidity Relative humidity [1]
     *  @returns Dewpoint temperature [K]
    */
    float dewpoint(float temperature, float relative_humidity);

    /** Vector version of dewpoint calculation
     *  @param temperature 1D vector of temperature [K]
     *  @param relative_humidity 1D vector of relative humidity [1]
     *  @returns 1D vector of dewpoint temperature [K]
    */
    vec dewpoint(const vec& temperature, const vec& relative_humidity);

    /** Calculate pressure at a new elevation
     *  @param ielev Elevation at start point [m]
     *  @param oelev Elevation at new point [m]
     *  @param ipressure Pressure at start point [Pa]
     *  @param itemperature Temperature at start point [K]
     *  @returns Pressure at new point [Pa]
     */
    float pressure(float ielev, float oelev, float ipressure, float itemperature=288.15);

    /** Calculate Vector version of pressure calculation
     *  @param ielev 1D vector of elevation at start point [m]
     *  @param oelev 1D vector of elevation at new point [m]
     *  @param ipressure 1D vector of pressure at start point [Pa]
     *  @param itemperature 1D vector of temperature at start point [K]
     *  @returns 1D vector of pressure at new points [Pa]
     */
    vec pressure(const vec& ielev, const vec& oelev, const vec& ipressure, const vec& itemperature);

    /** Convert Surface Pressure to Sea Level Pressure
     *  @param ps Surface pressure [Pa]
     *  @param altitude Station altitude above sea level [m]
     *  @param temperature 2m temperature [K]
     *  @param rh 2m Relative humidity [1]
     *  @param dewpoint 2m Dewpoint Temperature at station [K]
     *  @returns Sea Level Pressure [Pa]
     */
    float sea_level_pressure(float ps, float altitude, float temperature, float rh=gridpp::MV, float dewpoint=gridpp::MV);

    /** Vector version of Convert Surface Pressure to Sea Level Pressure
     *  @param ps 1D vector of Surface pressures [Pa]
     *  @param altitude 1D vector of station altitudes above sea level [m]
     *  @param temperature 1D vector of temperature [K]
     *  @param rh 1D vector of relative humidity [1]
     *  @param dewpoint 1D vector of dewpoint temperature [K]
     *  @returns 1D vector of sea Level pressure [Pa]
     */
    vec sea_level_pressure(const vec& ps, const vec& altitude, const vec& temperature, const vec& rh, const vec& dewpoint);

    /** Diagnose QNH from pressure and altitude
     *  @param pressure Pressure at point [Pa]
     *  @param altitude Altitude of point [m]
     *  @returns QNH [Pa]
    */
    float qnh(float pressure, float altitude);

    /** Vector version of QNH calculation
     *  @param pressure 1D vector of Pressure [Pa]
     *  @param altitude 1D vector of altitude [m]
     *  @returns 1D vector of QNH [Pa]
    */
    vec qnh(const vec& pressure, const vec& altitude);

    /** Calculate relative humidity from temperature and dewpoint temperature
     *  @param temperature Temperature [K]
     *  @param dewpoint Dewpoint temperature [K]
     *  @returns Relative humidity [1]
    */
    float relative_humidity(float temperature, float dewpoint);

    /** Vector version of relative humidity calculation
     *  @param temperature 1D vector of temperature [K]
     *  @param dewpoint 1D vector of dewpoint temperature [K]
     *  @returns 1D vector of relative humidity [1]
    */
    vec relative_humidity(const vec& temperature, const vec& dewpoint);

    /** Calculate wetbulb temperature from temperature, pressure, and relative humidity
     *  @param temperature Temperature [K]
     *  @param pressure Air pressure [Pa]
     *  @param relative_humidity Relative humidity [1]
     *  @returns Wetbulb temperature [K]
    */
    float wetbulb(float temperature, float pressure, float relative_humidity);

    /** Vector version of wetbulb calculation
     *  @param temperature 1D vector of temperature [K]
     *  @param pressure 1D vector of air pressure [pa]
     *  @param relative_humidity 1D vector of Relative humidity [1]
     *  @returns Wetbulb temperatures [K]
    */
    vec wetbulb(const vec& temperature, const vec& pressure, const vec& relative_humidity);

    /** Diagnose wind speed from its components
     *  @param xwind X-component of wind [any unit]
     *  @param ywind Y-component of wind [any unit]
     *  @returns Wind speed [any unit]
     * */
    float wind_speed(float xwind, float ywind);

    /** Vector version of wind speed calculation
     *  @param xwind 1D vector of X-component of wind [any unit]
     *  @param ywind 1D vector of Y-component of wind [any unit]
     *  @returns 1D vector of wind speeds [any unit]
     * */
    vec wind_speed(const vec& xwind, const vec& ywind);

    /** Diagnose wind direction from its components. If both xwind and ywind are 0, then direction
     *  is 180
     *  @param xwind X-component of wind [any unit]
     *  @param ywind Y-component of wind [any unit]
     *  @returns Wind direction [degrees]
     * */
    float wind_direction(float xwind, float ywind);

    /** Vector version of wind direction calculation
     *  @param xwind 1D vector of X-component of wind [any unit]
     *  @param ywind 1D vector of Y-component of wind [any unit]
     *  @returns 1D vector of wind direction [degrees]
     * */
    vec wind_direction(const vec& xwind, const vec& ywind);

    /**@}*/

    /** ****************************************
     * @name OpenMP settings
     * Functions that configure OpenMP
     * *****************************************/ /**@{*/
    /** Set the number of OpenMP threads to use. Overrides OMP_NUM_THREAD env variable.
     *  @param num Number of threads
    */
    void set_omp_threads(int num);

    /** Get the number of OpenMP threads currently set
     *  @returns Number of threads
    */
    int get_omp_threads();

    /** Sets the number of OpenMP threads to 1 if OMP_NUM_THREADS undefined */
    void initialize_omp();

    /** ****************************************
     * @name Utilities
     * Helper functions
     * *****************************************/ /**@{*/
    // ivec regression(const vec& x, const vec& y);

    /** Set the verbosity of debug messages. Use 0 for no messages.s
     *  @param level Debug level
    */
    void set_debug_level(int level);

    static int debug_level = 0;

    /** Get the currently set level of debug messagess
     *  @returns Currently set debug level
    */
    int get_debug_level();

    /** Convert name of a statistic enums
     *  @param name Name of a statistic
     *  @returns Statistic
    */
    Statistic get_statistic(std::string name);

    /** The gridpp version
     *  @returns The gridpp version
    */
    std::string version();

    /** The current time
     *  @returns The current time in seconds since 1870-01-01 00:00:00Z [s]
    */
    double clock();

    /** Writes a debug message to standard out
     *  @param string Write this string
    */
    void debug(std::string string);

    /** Writes a warning message to standard out
     *  @param string Write this string
    */
    void warning(std::string string);

    /** Writes an error message to standard out
     *  @param string Write this string
    */
    void error(std::string string);

    /** Writes an deprecation warning to standard out
     *  @param string Name of the deprecated function
     *  @param other Additional message to write
    */
    void future_deprecation_warning(std::string function, std::string other="");

    /** Checks if a value is valid
     *  @param value Value to check
     *  @returns True if the value is valid, False if it is NaN, Inf, or other invalid value
    */
    bool is_valid(float value);

    /** Compute a statistic on a 1D vector
     *  @param array 1D vector of values
     *  @param statistic Statistic to compute
     *  @returns Computed value
    */
    float calc_statistic(const vec& array, Statistic statistic);

    /** Compute a quantile from a 1D vector
     *  @param array 1D vector of values
     *  @param quantile_level Quantile level to compute
     *  @returns Computed quantile
    */
    float calc_quantile(const vec& array, float quantile_level);

    /** Compute a statistic on a 2D vector
     *  @param array 2D vector of values (Y, X)
     *  @param statistic Statistic to compute
     *  @returns 1D vector of computed values (Y)
    */
    vec calc_statistic(const vec2& array, Statistic statistic);

    /** Compute a quantile from a 2D vector
     *  @param array 2D vector of values (Y, X)
     *  @param quantile_level Quantile level to compute
     *  @returns 1D vector of computed quantiles (Y)
    */
    vec calc_quantile(const vec2& array, float quantile=MV);

    /** Compute quantile with 2D varying quantile level
      * @param array 3D input vector (T, Y, X)
      * @param quantile 2D vector of Quantile levels (Y, X)
      * @returns Extracted quantiles (Y, X)
    */
    vec2 calc_quantile(const vec3& array, const vec2& quantile_levels);

    /** Convert lats/lons or 2D x/y coordinates to 3D cartesian coordinates with the centre of the earth as the origin
     *  @param lats vector of latitudes [deg] or 2D y-coordinates [m]
     *  @param lons vector of longitudes [deg] or 2D x-coordinates [m]
     *  @param type Coordinate type for lats and lons
     *  @param x_coords vector of x-coordinates [m]
     *  @param y_coords vector of y-coordinates [m]
     *  @param z_coords vector of z-coordinates [m]
    */
    bool convert_coordinates(const vec& lats, const vec& lons, CoordinateType type, vec& x_coords, vec& y_coords, vec& z_coords);

    /** Convert lat/lon or 2D x/y coordinate to 3D cartesian coordinate with the centre of the earth as the origin
     *  @param lat vector of latitudes [deg] or 2D y-coordinates [m]
     *  @param lon vector of longitudes [deg] or 2D x-coordinates [m]
     *  @param type Coordinate type for lat and lon
     *  @param x_coord x-coordinate [m]
     *  @param y_coord y-coordinate [m]
     *  @param z_coord z-coordinate [m]
    */
    bool convert_coordinates(float lat, float lon, CoordinateType type, float& x_coord, float& y_coord, float& z_coord);

    /** Checks that a lat-coordinate is valid (based on the coordinate type)
     *  @param lat latitude [deg] or y-coordinate [m]
     *  @param type Coordinate type for lat
     *  @returns True if lat is valid latitude or y-coordinate
    */
    bool is_valid_lat(float lat, CoordinateType type);

    /** Checks that a lon-coordinate is valid (based on the coordinate type)
     *  @param lat longitude [deg] or x-coordinate [m]
     *  @param type Coordinate type for lon
     *  @returns True if lon is valid latitude or x-coordinate
    */
    bool is_valid_lon(float lon, CoordinateType type);

    /** Count the number of missing values in a vector
      * @param array 2D vector
      * @returns Number of invalid values
    */
    int num_missing_values(const vec2& iArray);

    /** Find the index in a vector that is equal to or just below a value
     *  @param iX Lookup value
     *  @param iValues 1D vector of lookup values. Must be sorted.
     *  @returns The index into iValues that is equal to or just below iX
    */
    int get_lower_index(float iX, const vec& iValues);

    /** Find the index in a vector that is equal to or just above a value
     *  @param iX Lookup value
     *  @param iValues 1D vector of lookup values. Must be sorted.
     *  @returns The index into iValues that is equal to or just above iX
    */
    int get_upper_index(float iX, const vec& iValues);

    /** Piecewise linear interpolation
     *  If x is outside the range of iX, then the min/max value of iY is used. If there are multiple
     *  identical x values, then the average of the y values at each end of the x-interval that is
     *  constant is used. The exception is if the constant interval is on one (only) of the edges of the
     *  interpolation curve. In that case, the y-value at the end of the interval further away from
     *  the boundary of the curve is used.
     *  @param x Interpolation to this point
     *  @param iX 1D vector of x-axis values. Vector must be sorted.
     *  @param iY 1D vector of y-axis values corresponding to iX.
     *  @returns Y value corresponding to x
    */
    float interpolate(float x, const vec& iX, const vec& iY);

    /** Piecewise linear interpolation, vectorised version.
     *  @param x 1D vector of values to interpolate to
     *  @param iX 1D vector of x-axis values. Vector must be sorted.
     *  @param iY 1D vector of y-axis values corresponding to iX.
     *  @returns 1D vector of y-axis values corresponding to x
    */
    vec interpolate(const vec& x, const vec& iX, const vec& iY);

    /** Initialize a 2D integer vector of size Y, X, with a given value
     *  @param Y Size of first dimension
     *  @param X Size of second dimension
     *  @patam value Initialize with this value
     *  @returns 2D interger vector
    */
    ivec2 init_ivec2(int Y, int X, int value);

    /** Initialize a 2D float vector of size Y, X, with a given value
     *  @param Y Size of first dimension
     *  @param X Size of second dimension
     *  @patam value Initialize with this value
     *  @returns 2D float vector
    */
    vec2 init_vec2(int Y, int X, float value=MV);

    /** Initialize a 3D integer vector of size Y, X, E, with a given value
     *  @param Y Size of first dimension
     *  @param X Size of second dimension
     *  @param E Size of third dimension
     *  @patam value Initialize with this value
     *  @returns 3D integer vector
    */
    ivec3 init_ivec3(int Y, int X, int E, int value);

    /** Initialize a 3D float vector of size Y, X, E, with a given value
     *  @param Y Size of first dimension
     *  @param X Size of second dimension
     *  @param E Size of third dimension
     *  @patam value Initialize with this value
     *  @returns 3D float vector
    */
    vec3 init_vec3(int Y, int X, int E, float value=MV);

    /** Get reasonably spaced quantile levels from a vector of values, ignoring duplicate values
      *  but including the first number after duplicated values. Include the lowest and highest
      *  values.
      *  @param values 1D vector of values (unsorted, and no invalid values)
      *  @param num number of thresholds to get
      *  @returns 1D vector of sorted quantile levels
    */
    vec calc_even_quantiles(const vec& values, int num);

    /** Compute window statistics across time, independently for each case
    *  @param array 2D vector of input values (Case, Time)
    *  @param length Window length in number of timesteps
    *  @param statistic Statistic to apply to window
    *  @param before If true, make the window end at the particular time. If false, centre it.
    *  @param keep_missing If true, window value will be missing if one or more values in window are missing
    *  @param missing_edges If true put missing values at the edges, where window overshoots the edge
    *  @returns 2D vector of output values (Case, Time)
    */
    vec2 window(const vec2& array, int length, gridpp::Statistic statistic, bool before=false, bool keep_missing=false, bool missing_edges=true);

    /** Check if the grid is the same size as the 2D vector
     *  @param grid Grid (Y, X)
     *  @param v 2D vector (Y, X)
     *  @returns True if compatible
    */
    bool compatible_size(const Grid& grid, const vec2& v);

    /** Check if the grid is compatible with 3D vector
     *  @param grid Grid (Y, X)
     *  @param v 3D vector (T, Y, X)
     *  @returns True if compatible
    */
    bool compatible_size(const Grid& grid, const vec3& v);

    /** Check if the points are compatible with a 1D vector
     *  @param grid Points
     *  @param v 1D vector
     *  @returns True if compatible
    */
    bool compatible_size(const Points& points, const vec& v);

    /** Check if the points are compatible with a 2D vector
     *  @param grid Points (P)
     *  @param v 2D vector (T, P)
     *  @returns True if compatible
    */
    bool compatible_size(const Points& points, const vec2& v);

    /** Check if two 2D vectors are compatible
     *  @param a 2D vector (Y, X)
     *  @param b 2D vector (Y, X)
     *  @returns True if compatible
    */
    bool compatible_size(const vec2& a, const vec2& b);

    /** Check if a 2D and 3D vectors are compatible
     *  @param a 2D vector (Y, X)
     *  @param b 3D vector (Y, X, E)
     *  @returns True if compatible
    */
    bool compatible_size(const vec2& a, const vec3& b);

    /** Check if two 3D vectors are compatible
     *  @param a 3D vector (Y, X, E)
     *  @param b 3D vector (Y, X, E)
     *  @returns True if compatible
    */
    bool compatible_size(const vec3& a, const vec3& b);

    /** Checks if a point is located inside a rectangle formed by 4 points. The 4 points must be
      * provided in an order that draws out a rectangle (either clockwise or counter-clockwise).
      * @param A: A point in the rectangle
      * @param B: A point in the rectangle
      * @param C: A point in the rectangle
      * @param D: A point in the rectangle
      * @param m: The point to test if it is inside
      * @returns True if the point is inside, False otherwise
    */
    bool point_in_rectangle(const Point& A, const Point& B, const Point& C, const Point& D, const Point& m );

    /** ****************************************
     * @name SWIG testing functions
     * Functions for testing the SWIG interface. Not useful for any other purpose.
     * *****************************************/ /**@{*/
    /** Special function whose presense is needed for SWIG */
    float* test_array(float* v, int n);
    /** Testing function for 1D input vector */
    float test_vec_input(const vec& input);
    /** Testing function for 1D input vector */
    int test_ivec_input(const ivec& input);
    /** Testing function for 2D input vector */
    float test_vec2_input(const vec2& input);
    /** Testing function for 3D input vector */
    float test_vec3_input(const vec3& input);
    /** Testing function for 1D output vector */
    vec test_vec_output();
    ivec test_ivec_output();
    /** Testing function for 2D output vector */
    vec2 test_vec2_output();
    ivec2 test_ivec2_output();
    /** Testing function for 3D output vector */
    vec3 test_vec3_output();
    ivec3 test_ivec3_output();

    /** Testing function for 1D vector treated as output */
    float test_vec_argout(vec& distances);
    /** Testing function for 2D vector treated as output */
    float test_vec2_argout(vec2& distances);

    void test_not_implemented_exception();
/*    vec2 test_args_for_R(const Points& bpoints, const StructureFunction& structure, const vec2& background);*/
    vec2 test_args_for_R(const Points& bpoints, const vec2& background);

    /** Default value used to fill array in SWIG testing functions. Not useful for any other purpose. */
    static const float swig_default_value = -1;

    /**@}*/

/** Represents a single point in some coordinate system */
    class Point {
        public:
            /** Constructor
              * @param lat: Latitude (or y) coordinate
              * @param lon: Longitude (or x) coordinate
              * @param elev: Elevation
              * @param laf: Land area fraction (between 0 and 1)
              * @param type: Coordinate type for lat and lon
            */
            Point(float lat, float lon, float elev=MV, float laf=MV, CoordinateType type=Geodetic);

            /** Constructor
              * @param lat: Latitude (or y) coordinate
              * @param lon: Longitude (or x) coordinate
              * @param elev: Elevation
              * @param laf: Land area fraction (between 0 and 1)
              * @param type: Coordinate type for lat and lon
              * @param x: X-coordinate
              * @param y: Y-coordinate
              * @param x: Z-coordinate
            */
            Point(float lat, float lon, float elev, float laf, CoordinateType type, float x, float y, float z);
            float lat;
            float lon;
            float elev;
            float laf;
            CoordinateType type;
            float x;
            float y;
            float z;
    };

    /** Helper class for Grid and Points representing a tree of points */
    class KDTree {
        public:
            KDTree(vec lats, vec lons, CoordinateType type=Geodetic);
            KDTree& operator=(KDTree other);
            KDTree(const KDTree& other);
            KDTree(CoordinateType type=Geodetic) {mType=type;};

            /** Find single nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             * */
            int get_nearest_neighbour(float lat, float lon, bool include_match=true) const;

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            ivec get_neighbours(float lat, float lon, float radius, bool include_match=true) const;

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             * */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match=true) const;

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            int get_num_neighbours(float lat, float lon, float radius, bool include_match=true) const;

            /** Find a set of nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param num Number of points to find
             * */
            ivec get_closest_neighbours(float lat, float lon, int num, bool include_match=true) const;

            static float deg2rad(float deg);
            static float rad2deg(float deg);

            /** Calculate distance between two coordinates using Haversine function
             *  @param lat1 Latitude [deg] or y-coordinate [m] for first point
             *  @param lon1 Longitude [deg] or x-coordinate [m] for first point
             *  @param lat2 Latitude [deg] or y-coordinate [m] for second point
             *  @param lon2 Longitude [deg] or x-coordinate [m] for second point
             *  @param type Coordinate type for both points
             *  @returns Distance [m]
            **/
            static float calc_distance(float lat1, float lon1, float lat2, float lon2, CoordinateType type=Geodetic);

            /** Calculate straight line distance between two 3D coordinates. This is much faster
             *  than calc_distance.
             *  @param x0 x-coordinate of first point [m]
             *  @param y0 y-coordinate of first point [m]
             *  @param z0 z-coordinate of first point [m]
             *  @param x1 x-coordinate of second point [m]
             *  @param y1 y-coordinate of second point [m]
             *  @param z1 z-coordinate of second point [m]
             *  @returns Distance [m]
            **/
            static float calc_straight_distance(float x0, float y0, float z0, float x1, float y1, float z1);

            /** Calculate distance between two points using Haversine function
             *  @param p1 First point
             *  @param p2 Second point
             *  @returns Distance [m]
            **/
            static float calc_distance(const Point& p1, const Point& p2);

            /** Calculate straight line distance between two points. This is much faster than
             *  calc_distance.
             *  @param p1 First point
             *  @param p2 Second point
             *  @returns Distance [m]
            **/
            static float calc_straight_distance(const Point& p1, const Point& p2);

            /** Calculate an approximate distance between two coordinates
             *  @param lat1 Latitude [deg] or y-coordinate [m] for first point
             *  @param lon1 Longitude [deg] or x-coordinate [m] for first point
             *  @param lat2 Latitude [deg] or y-coordinate [m] for second point
             *  @param lon2 Longitude [deg] or x-coordinate [m] for second point
             *  @param type Coordinate type for both points
             *  @returns Distance [m]
            **/
            static float calc_distance_fast(float lat1, float lon1, float lat2, float lon2, CoordinateType type=Geodetic);

            vec get_lats() const;
            vec get_lons() const;
            int size() const;
            CoordinateType get_coordinate_type() const;
            const vec& get_x() const;
            const vec& get_y() const;
            const vec& get_z() const;
        protected:
            typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
            typedef std::pair<point, unsigned> value;
            typedef boost::geometry::model::box<point> box;
            boost::geometry::index::rtree< value, boost::geometry::index::quadratic<16> > mTree;
            vec mLats;
            vec mLons;
            vec mX;
            vec mY;
            vec mZ;
            CoordinateType mType;

            struct within_radius {
                public:
                    within_radius(point p, float radius, bool include_match);
                    bool operator()(value const& v) const;
                private:
                    float radius;
                    point p;
                    bool include_match;
            };
            struct is_not_equal {
                public:
                    is_not_equal(point p);
                    bool operator()(value const& v) const;
                private:
                    point p;
            };
    };

    /** Represents a vector of locations and their metadata */
    class Points  {
        public:
            Points();
            /** Initialize a new grid
             *  @param lats: 1D vector of latitudes [degrees]
             *  @param lons: 1D vector of longitudes [degrees]
             *  @param elevs: 1D vector of elevations [m]
             *  @param lafs: 1D vector of land area fractions [1]
             *  @param type: Coordinate type
            */
            Points(vec lats, vec lons, vec elevs=vec(), vec lafs=vec(), CoordinateType type=Geodetic);
            Points(KDTree tree, vec elevs=vec(), vec lafs=vec());
            Points& operator=(Points other);
            Points(const Points& other);

            /* Get the index of the nearest neighbour
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param include_match Return matching points. If false, don't use this point.
             * @returns Index of nearest neighbour (-1 if there are no neighbours).
            */
            int get_nearest_neighbour(float lat, float lon, bool include_match=true) const;

            /* Get the indices of all neighbours within a search radius
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param radius Neighbourhood radius [m]
             * @param include_match Return matching points. If false, don't use this point.
             * @returns 1D vector of indices of neighbours within radius
            */
            ivec get_neighbours(float lat, float lon, float radius, bool include_match=true) const;

            /* Get the indices and distances to all neighbours within a search radius
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param radius Neighbourhood radius [m]
             * @param distances 1D output vector of distances to each point [m]
             * @param include_match Return matching points. If false, don't use this point.
             * @returns 1D vector of indices of neighbours within radius
            */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match=true) const;

            /* Get the number of neighbouring points within a search radius
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param radius Neighbourhood radius [m]
             * @param include_match Return matching points. If false, don't use this point.
             * @returns the number of neighbours within radius
            */
            int get_num_neighbours(float lat, float lon, float radius, bool include_match=true) const;

            /* Get the N nearest neighbouring points
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param num Number of nearest neighbours
             * @param include_match Return matching points. If false, don't use this point.
             * @returns 1D vector of indices of neighbours
            */
            ivec get_closest_neighbours(float lat, float lon, int num, bool include_match=true) const;

            vec get_lats() const;
            vec get_lons() const;
            vec get_elevs() const;
            vec get_lafs() const;

            int size() const;

            /** Get point indices that are inside the envelope of a grid
             *  @param grid Grid
             *  @returns 1D vector of indices
            */
            ivec get_in_domain_indices(const Grid& grid) const;

            /** Get points that are inside the envelope of a grid
             *  @param grid Grid
             *  @returns A new points object containing points within grid
            */
            Points get_in_domain(const Grid& grid) const;
            CoordinateType get_coordinate_type() const;
            Point get_point(int index) const;

            /** Subset the points
             *  @param indices 1D vector of point indices
             *  @returns A new points object containing the indexed points
            */
            Points subset(const ivec& indices) const;
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

            /** Initialize a new grid
             *  @param lats: 2D vector of latitudes [degrees] (or y-values)
             *  @param lons: 2D vector of longitudes [degrees] (or x-values)
             *  @param elevs: 2D vector of elevations [m]
             *  @param lafs: 2D vector of land area fractions [1]
             *  @param type: Coordinate type
            */
            Grid(vec2 lats, vec2 lons, vec2 elevs=vec2(), vec2 lafs=vec2(), CoordinateType type=Geodetic);

            /* Get the index of the nearest neighbour
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param include_match Return matching points. If false, don't use this point.
             * @returns 1D vector of (x, y) index (empty vector if no neighbours)
            */
            ivec get_nearest_neighbour(float lat, float lon, bool include_match=true) const;

            /* Get the indices of all neighbours within a search radius
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param radius Neighbourhood radius [m]
             * @param include_match Return matching points. If false, don't use this point.
             * @returns 2D vector (Point, 2) of x,y-indices of neighbours within radius
            */
            ivec2 get_neighbours(float lat, float lon, float radius, bool include_match=true) const;

            /* Get the indices and distances to all neighbours within a search radius
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param radius Neighbourhood radius [m]
             * @param distances 1D output vector of distances to each point [m]
             * @param include_match Return matching points. If false, don't use this point.
             * @returns 2D vector (Point, 2) of x,y-indices of neighbours within radius
            */
            ivec2 get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match=true) const;

            /* Get the number of neighbouring points within a search radius
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param radius Neighbourhood radius [m]
             * @param include_match Return matching points. If false, don't use this point.
             * @returns the number of neighbours within radius
            */
            int get_num_neighbours(float lat, float lon, float radius, bool include_match=true) const;

            /* Get the N nearest neighbouring points
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param num Number of nearest neighbours
             * @param include_match Return matching points. If false, don't use this point.
             * @returns 2D vector (Point, 2) of x,y-indices of neighbours
            */
            ivec2 get_closest_neighbours(float lat, float lon, int num, bool include_match=true) const;

            /* Get the N nearest neighbouring points
             * @param lat Lookup latitude (or y-axis value)
             * @param lon Lookup longitude (or x-axis value)
             * @param Y1_out Number of nearest neighbours
             * @param include_match Return matching points. If false, don't use this point.
             * @returns 2D vector (Point, 2) of x,y-indices of neighbours
            */
            bool get_box(float lat, float lon, int& Y1_out, int& X1_out, int& Y2_out, int& X2_out) const;

            /** Convert grid to a vector of points
             *  @returns Points object containing all gridpoints serialised into a vector
            */
            Points to_points() const;

            vec2 get_lats() const;
            vec2 get_lons() const;
            vec2 get_elevs() const;
            vec2 get_lafs() const;
            ivec size() const;
            CoordinateType get_coordinate_type() const;
            Point get_point(int y_index, int x_index) const;
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
    class not_implemented_exception: public std::logic_error
    {
        public:
            not_implemented_exception();
    };

    typedef std::shared_ptr<gridpp::StructureFunction> StructureFunctionPtr;
    /** Covariance structure function */
    class StructureFunction {
        public:
            StructureFunction(float localization_distance=0);
            /** Correlation between two points
             *  @param p1 First point
             *  @param p2 Other point
             *  @returns Correlation between points
            */
            virtual float corr(const Point& p1, const Point& p2) const = 0;
            virtual vec corr(const Point& p1, const std::vector<Point>& p2) const;

            /** Correlation between a background point and an observation points
             *  @param p1 Background point
             *  @param p2 Observation point
             *  @returns Correlation between background and observation points
            */
            virtual float corr_background(const Point& p1, const Point& p2) const;
            virtual vec corr_background(const Point& p1, const std::vector<Point>& p2) const;

            /** Maximum distance for which an observation can have an impact (localization)
              * @returns Distance [m]
            */
            virtual float localization_distance(const Point& p) const;
            virtual StructureFunctionPtr clone() const = 0;
            static const float default_min_rho;
        protected:
            /** Barnes correlation function
              * @param dist Distance between points. Same units as 'length'
              * @param length Length scale
              * @returns Barnes rho
            */
            float barnes_rho(float dist, float length) const;

            /** Cressman correlation function
              * @param dist Distance between points. Same units as 'length'
              * @param length Length scale
              * @returns Cressman rho
            */
            float cressman_rho(float dist, float length) const;

            /** Compactly supported second-order autoregressive correlation function
              * @param dist Distance between points. Same units as 'length'
              * @param length Length scale
              * @returns SOAR rho
            */
            float soar_rho(float dist, float length) const;

            /** Compactly supported third-order autoregressive correlation function
              * @param dist Distance between points. Same units as 'length'
              * @param length Length scale
              * @returns TOAR rho
            */
            float toar_rho(float dist, float length) const;

            /** Powerlaw correlation function
              * @param dist Distance between points. Same units as 'length'
              * @param length Length scale
              * @returns powerlaw rho
            */
            float powerlaw_rho(float dist, float length) const;

            /** Linear correlation function
              * @param normdist Normalized distance between points. Must be in the range -1 to 1.
              * @param min Minimum allowed value for the correlation (if less than 0, the return 1)
              * @returns linear rho
            */
            float linear_rho(float normdist, float min) const;
            float m_localization_distance;
    };
    class MultipleStructure: public StructureFunction {
        public:
            /** Different structure functions for horizontal, vertical, and land/sea
              * @param structure_h: Horizontal structure function
              * @param structure_v: Vertical structure function
              * @param structure_w: Land/sea structure function
            */
            MultipleStructure(const StructureFunction& structure_h, const StructureFunction& structure_v, const StructureFunction& structure_w);
            float corr(const Point& p1, const Point& p2) const;
            vec corr(const Point& p1, const std::vector<Point>& p2) const;
            StructureFunctionPtr clone() const;
            float localization_distance(const Point& p) const;
        private:
            StructureFunctionPtr m_structure_h;
            StructureFunctionPtr m_structure_v;
            StructureFunctionPtr m_structure_w;
    };
    /** Simple structure function based on distance, elevation, and land area fraction */
    class BarnesStructure: public StructureFunction {
        public:
            /** Exponential structure function
              * @param h: Horizontal decorrelation length >=0 [m]
              * @param v: Vertical decorrelation length >=0 [m]. If 0, disable decorrelation.
              * @param w: Land/sea decorrelation length >=0 [1]. If 0, disable decorrelation.
              * @param hmax: Truncate horizontal correlation beyond this length [m]. If undefined, 3.64 * h.
            */
            BarnesStructure(float h, float v=0, float w=0, float hmax=MV);

            /** Barnes structure function where decorrelation varyies spatially
              * @param grid: Grid of decorrelation field
              * @param h: 2D vector of horizontal decorrelation lengths >=0, same size as grid [m]
              * @param v: 2D vector of Vertical decorrelation lengths >=0 [m]. Set all to 0 to disable decorrelation.
              * @param w: 2D vector of land/sea decorrelation lengths >=0 [1]. Set all to 0 to disable decorrelation.
              * @param min_rho: Truncate horizontal correlation when rho is less than this value [m].
            */
            BarnesStructure(Grid grid, vec2 h, vec2 v, vec2 w, float min_rho=StructureFunction::default_min_rho);
            float corr(const Point& p1, const Point& p2) const;
            vec corr(const Point& p1, const std::vector<Point>& p2) const;
            StructureFunctionPtr clone() const;
            float localization_distance(const Point& p) const;
        private:
            float localization_distance(float h) const;
            Grid m_grid;
            vec2 mH;
            vec2 mV;
            vec2 mW;
            float m_min_rho;
            bool m_is_spatial;
    };

    /** Mix of structure functions based on distance(SOAR), elevation(Gauss), and land area fraction(linear) */
    class MixAStructure: public StructureFunction {
        public:
            /** Structure function SOAR - Gaussian - Linear
              * @param h: Horizontal decorrelation length >=0 [m] (SOAR)
              * @param v: Vertical decorrelation length >=0 [m] (Gaussian). If 0, disable decorrelation.
              * @param w: Land/sea decorrelation length >=0 [1] (Linear). If 0, disable decorrelation.
              * @param hmax: Truncate horizontal correlation beyond this length [m]. If undefined, 3.64 * h.
            */
            MixAStructure(float h, float v=0, float w=0, float hmax=MV);

            /** MixA structure function where decorrelation varyies spatially
              * @param grid: Grid of decorrelation field
              * @param h: 2D vector of horizontal decorrelation lengths >=0 (SOAR) , same size as grid [m]
              * @param v: 2D vector of Vertical decorrelation lengths >=0 [m] (Gaussian). Set all to 0 to disable decorrelation.
              * @param w: 2D vector of land/sea decorrelation lengths >=0 [1] (Linear). Set all to 0 to disable decorrelation.
              * @param min_rho: Truncate horizontal correlation when rho is less than this value [m].
            */
            MixAStructure(Grid grid, vec2 h, vec2 v, vec2 w, float min_rho=StructureFunction::default_min_rho);
            float corr(const Point& p1, const Point& p2) const;
            vec corr(const Point& p1, const std::vector<Point>& p2) const;
            StructureFunctionPtr clone() const;
            float localization_distance(const Point& p) const;
        private:
            float localization_distance(float h) const;
            Grid m_grid;
            vec2 mH;
            vec2 mV;
            vec2 mW;
            float m_min_rho;
            bool m_is_spatial;
    };

    /** Simple structure function based on distance, elevation, and land area fraction */
    class CressmanStructure: public StructureFunction {
        public:
            CressmanStructure(float h, float v=0, float w=0);
            float corr(const Point& p1, const Point& p2) const;
            StructureFunctionPtr clone() const;
        private:
            float mH;
            float mV;
            float mW;
    };

    class CrossValidation: public StructureFunction {
        public:
            /** Structure function for performing cross validation experiments
              * @param structure: Base structure function to use
              * @param dist: Force background-to-obs correlation to 0 for points within this distance [m]
            */
            CrossValidation(StructureFunction& structure, float dist);
            float corr(const Point& p1, const Point& p2) const;
            float corr_background(const Point& p1, const Point& p2) const;
            vec corr_background(const Point& p1, const std::vector<Point>& p2) const;
            StructureFunctionPtr clone() const;
            float localization_distance(const Point& p) const;
        private:
            StructureFunctionPtr m_structure;
            float m_dist;
    };

    class Transform {
        public:
            // Note these cannot be pure virtual, otherwise SWIG does not expose
            // the vector functions (with the same name) in python. Therefore, make sure these
            // functions are overloaded in the subclass implementation
            virtual float forward(float value) const;
            virtual float backward(float value) const;

            /** Forward transform a 1D vector
             *  @param input 1D vector
             *  @returns 1D vector of transformed values
            */
            vec forward(const vec& input) const;

            /** Backward transform a 1D vector
             *  @param input 1D vector of transformed values
             *  @returns 1D vector of original values
            */
            vec backward(const vec& input) const;

            /** Forward transform a 2D vector
             *  @param input 2D vector
             *  @returns 2D vector of transformed values
            */
            vec2 forward(const vec2& input) const;

            /** Backward transform a 2D vector
             *  @param input 2D vector of transformed values
             *  @returns 2D vector of original values
            */
            vec2 backward(const vec2& input) const;

            /** Forward transform a 3D vector
             *  @param input 3D vector
             *  @returns 3D vector of transformed values
            */
            vec3 forward(const vec3& input) const;

            /** Backward transform a 3D vector
             *  @param input 3D vector of transformed values
             *  @returns 3D vector of original values
            */
            vec3 backward(const vec3& input) const;
    };
    /** Identity transform, i.e. where forward and backward functinos to not modify values */
    class Identity : public Transform {
        public:
            // SWIG requires these "using" statements to enable the vectorized versions in the
            // subclasses
            using Transform::forward;
            using Transform::backward;
            float forward(float value) const;
            float backward(float value) const;
    };
    /** Log transformation: output = log(input) */
    class Log : public Transform {
        public:
            using Transform::forward;
            using Transform::backward;
            float forward(float value) const;
            float backward(float value) const;
    };
    /** Box-Cox transformation */
    class BoxCox : public Transform {
        public:
            /** Initialize Box-Cox transform
             *  @param threshold Box-Cox parameter
            */
            BoxCox(float threshold);
            using Transform::forward;
            using Transform::backward;
            float forward(float value) const;
            float backward(float value) const;
        private:
            float mThreshold;
    };
    /** Gamma transformation. Transforms values to cdf from a gamma distribution and subsequantly
     *  extracts the cdf from a standard normal distribution. */
    class Gamma : public Transform {
        public:
            /** Initialize Gamma transform
             *  @param shape Shape parameter
             *  @param scale Scale parameter
             *  @param tolerance Computation tolerance
             *
            */
            Gamma(float shape, float scale, float tolerance=0.01);
            using Transform::forward;
            using Transform::backward;
            float forward(float value) const;
            float backward(float value) const;
        private:
            float m_tolerance;
            boost::math::gamma_distribution<> m_gamma_dist;
            boost::math::normal m_norm_dist;
    };
};
#endif
