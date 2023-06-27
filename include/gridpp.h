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

#define GRIDPP_VERSION "0.7.0.dev4"
#define __version__ GRIDPP_VERSION

namespace gridpp {
    /** **************************************
     * @name Short-hand notation for vectors of different dimensions sizes
     * ***************************************/ /**@{*/
    // Preferred vector types
    typedef std::vector<float> vec;
    typedef std::vector<vec> vec2;
    typedef std::vector<vec2> vec3;
    typedef std::vector<int> ivec;
    typedef std::vector<ivec> ivec2;
    typedef std::vector<ivec2> ivec3;

    // Only use when double is required
    typedef std::vector<double> dvec;
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
        RandomChoice = 90, /**< Randomly pick a value */
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
        Nearest = 0,
        Bilinear = 1,
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
      * @param points Points of observations
      * @param pobs Vector of observations
      * @param pratios Vector of ratio of observation error variance to background variance
      * @param pbackground Background with observation operator
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
    */
    vec2 optimal_interpolation(const Grid& bgrid,
            const vec2& background,
            const Points& points,
            const vec& pobs,
            const vec& pratios,
            const vec& pbackground,
            const StructureFunction& structure,
            int max_points,
            bool allow_extrapolation=true);

    /** Optimal interpolation for a deterministic vector of points
      * @param bpoints Points of background field
      * @param background 1D field of background values
      * @param points Points of observations
      * @param pobs Vector of observations
      * @param pratios Vector of ratio of observation error variance to background variance
      * @param pbackground Background with observation operator
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
    */
    vec optimal_interpolation(const Points& bpoints,
            const vec& background,
            const Points& points,
            const vec& pobs,
            const vec& pratios,
            const vec& pbackground,
            const StructureFunction& structure,
            int max_points,
            bool allow_extrapolation=true);

    /** Optimal interpolation for a deterministic gridded field including analysis variance
      * @param bgrid Grid of background field
      * @param background 2D field of background values
      * @param bvariance Variance of background field
      * @param points Points of observations
      * @param obs Vector of observations
      * @param obs_variance Variance of observations
      * @param background_at_points Background interpolated to observation points
      * @param bvariance_at_points Variance of background interpolated to observation points
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
    */
    vec2 optimal_interpolation_full(const Grid& bgrid,
            const vec2& background,
            const vec2& bvariance,
            const Points& points,
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
      * @param background 1D field of background values
      * @param bvariance Variance of background field
      * @param points Points of observations
      * @param obs Vector of observations
      * @param obs_variance Variance of observations
      * @param background_at_points Background interpolated to observation points
      * @param bvariance_at_points Variance of background interpolated to observation points
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param allow_extrapolation Allow OI to extrapolate increments outside increments at observations
    */
    vec optimal_interpolation_full(const Points& bpoints,
            const vec& background,
            const vec& bvariance,
            const Points& points,
            const vec& obs,
            const vec& obs_variance,
            const vec& background_at_points,
            const vec& bvariance_at_points,
            const StructureFunction& structure,
            int max_points,
            vec& analysis_sigmas,
            bool allow_extrapolation=true);

    /** Optimal interpolation using a structure function based on an ensemble 
      * See Lussana et al 2019 (DOI: 10.1002/qj.3646)
      * @param input 3D field of background values (Y, X, E)
      * @param bgrid grid corresponding to input
      * @param pobs vector of observations
      * @param pci vector of ci values
      * @param points observation points
    */
    vec3 optimal_interpolation_ensi(const Grid& bgrid,
            const vec3& background,
            const Points& points,
            const vec& pobs,
            const vec& psigmas,
            const vec2& pbackground,
            const StructureFunction& structure,
            int max_points,
            bool allow_extrapolation=true);

    vec2 optimal_interpolation_ensi(const Points& bpoints,
            const vec2& background,
            const Points& points,
            const vec& pobs,
            const vec& psigmas,
            const vec2& pbackground,
            const StructureFunction& structure,
            int max_points,
            bool allow_extrapolation=true);

    /** Correction of a gridded field ensuring the distribution of values nearby match that of
      * observations. This is an experimental method.
      * @param bgrid grid corresponding to input
      * @param background 2D field of background values (Y, X)
      * @param points observation points
      * @param pobs vector of observations
      * @param pbackground vector of background values at points
      * @param structure structure function specifying correlation between points
      * @param min_quantile truncate quantile map below this quantile
      * @param max_quantile truncate quantile map above this quantile
      * @param max_points maximum number of points used within localization radius (not used at moment)
    */
    vec2 local_distribution_correction(const Grid& bgrid,
            const vec2& background,
            const Points& points,
            const vec& pobs,
            const vec& pbackground,
            const StructureFunction& structure,
            float min_quantile,
            float max_quantile,
            int min_points=0);

    /** Version with multiple number of timesteps. Pool all observations across time in when
      * computing the calibration curve.
      * @param bgrid grid corresponding to input
      * @param background 2D field of background values (Y, X)
      * @param points observation points
      * @param pobs 2D vector of observations with dimensions (T, N)
      * @param pbackground vector of background values at points with dimensions (T, N)
      * @param structure structure function specifying correlation between points
      * @param min_quantile truncate quantile map below this quantile
      * @param max_quantile truncate quantile map above this quantile
      * @param max_points maximum number of points used within localization radius (not used at moment)
    */
    vec2 local_distribution_correction(const Grid& bgrid,
            const vec2& background,
            const Points& points,
            const vec2& pobs,
            const vec2& pbackground,
            const StructureFunction& structure,
            float min_quantile,
            float max_quantile,
            int min_points=0);

    /** Fill in values inside or outside a set of circles (useful for masking)
      * @param input Deterministic values with dimensions Y, X
      * @param radii Circle radii for each point
      * @param value Fill in this value
      * @param outside if True, fill outside circles, if False, fill inside circles
    */
    vec2 fill(const Grid& igrid, const vec2& input, const Points& points, const vec& radii, float value, bool outside);


    /** Insert observations into gridded field using a square box
      * @param grid Grid
      * @param background Deterministic values with dimensions Y, X
      * @param points Points representing observations
      * @param observations Vector of observations
      * @param halfwidths Half width of square (in number of grid points) where observations are inserted for each point
      * @param max_elev_diff Only insert where elevation difference between grid and point is less than this value
    */
    vec2 doping_square(const Grid& igrid, const vec2& background, const Points& points, const vec& observations, const ivec& halfwidths, float max_elev_diff=gridpp::MV);

    /** Insert observations into gridded field using a circle
      * @param grid Grid
      * @param background Deterministic values with dimensions Y, X
      * @param points Points representing observations
      * @param observations Vector of observations
      * @param radii Radius of circle where observations are inserted for each point [m]
      * @param max_elev_diff Only insert where elevation difference between grid and point is less than this value
    */
    vec2 doping_circle(const Grid& igrid, const vec2& background, const Points& points, const vec& observations, const vec& radii, float max_elev_diff=gridpp::MV);

    /** **************************************
     * @name Distributions
     * Functions that extract values from probability distributions
     * ***************************************/ /**@{*/

    /**
     * @param shape Shape parameter of gamma distribution
     * @param scale Scale parameter of gamma distribution
     * @param levels Quantile levels to retrieve
    */
    vec gamma_inv(const vec& levels, const vec& shape, const vec& scale);

    /**@}*/

    /** **************************************
     * @name Spatial neighbourhood filters
     * Functions that apply neighbourhood filters on a gridded field
     * ***************************************/ /**@{*/

    /** Spatial neighbourhood filter, computing a statistic for a sliding square window
      * @param input 2D grid of values
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param statistic Statistic to compute
    */
    vec2 neighbourhood(const vec2& input, int halfwidth, Statistic statistic);

    /** Spatial neighbourhood filter for an ensemble of fields
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param statistic Statistic to compute
    */
    vec2 neighbourhood(const vec3& input, int halfwidth, Statistic statistic);

    /** Computes a quantile in a sliding square neighbourhood
      * @param input 2D grid of values
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
    */
    vec2 neighbourhood_quantile(const vec2& input, float quantile, int halfwidth);

    /** Computes a quantile in a sliding square neighbourhood for an ensemble of fields
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
    */
    vec2 neighbourhood_quantile(const vec3& input, float quantile, int halfwidth);

    /** Fast and approximate neighbourhood quantile
      * @param input 2D grid of values
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec2& input, float quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile for ensemble of fields
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec3& input, float quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile, with spatially varying quantile
      * @param input 2D grid of values
      * @param quantile 2D grid quantiles to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec2& input, const vec2& quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile for ensemble of fields, with spatially varying quantile
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param quantile 2D grid quantiles to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec3& input, const vec2& quantile, int halfwidth, const vec& thresholds);

    /** Spatial neighbourhood filter without any shortcuts. This is quite slow and is only useful for testing.
      * @param input 2D grid of values
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param operation one of min, mean, median, max
    */
    vec2 neighbourhood_brute_force(const vec2& input, int halfwidth, Statistic statistic);

    /** Spatial neighbourhood filter without any shortcuts. This is quite slow and is only useful for testing.
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param operation one of min, mean, median, max
    */
    vec2 neighbourhood_brute_force(const vec3& input, int halfwidth, Statistic statistic);

    /** Calculate appropriate approximation thresholds for neighbourhood quantile
      * @param input 2D (Y, X) array of values
      * @param num_thresholds Number of thresholds
    */
    vec get_neighbourhood_thresholds(const vec2& input, int num_thresholds);

    /** Calculate appropriate approximation thresholds for neighbourhood quantile based on an * ensemble
      * @param input 3D (Y, X, T) array of values
      * @param num_thresholds Number of thresholds
    */
    vec get_neighbourhood_thresholds(const vec3& input, int num_thresholds);

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
     *  @param ref Reference values (observations)
     *  @param fcst Forecast values
     *  @param output_fcst Output forecast quantiles
     *  @param quantiles Vector of quantiles to extract. If empty, use all values.
     *  @return Output reference quantiles
    */
    vec quantile_mapping_curve(const vec& ref, const vec& fcst, vec& output_fcst, vec quantiles=vec());

    /** Create calibration curve that optimizes a metric
     *  @param ref Reference values (observations)
     *  @param fcst Forecast values
     *  @param thresholds Thresholds for computing optimal values for
     *  @param metric A Metric to optimize for
     *  @return Calibration curve
    */

    vec metric_optimizer_curve(const vec& ref, const vec& fcst, const vec& thresholds, Metric metric, vec& output_fcst);

    /** Apply arbitrary calibration curve to a single value
     *  @param fcst input forecast
     *  @param curve_ref Reference quantiles
     *  @param curve_fcst Forecast quantiles
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @return Calibrated forecasts
    */
    float apply_curve(float fcst, const vec& curve_ref, const vec& curve_fcst, Extrapolation policy_below, Extrapolation policy_above);

    /** Apply arbitrary calibration curve to 1D forecasts
     *  @param fcst 1D vector of forecast values
     *  @param curve_ref Reference quantiles
     *  @param curve_fcst Forecast quantiles
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @return Calibrated forecasts
    */
    vec apply_curve(const vec& fcst, const vec& curve_ref, const vec& curve_fcst, Extrapolation policy_below, Extrapolation policy_above);

    /** Apply arbitrary calibration curve to 2D forecasts
     *  @param fcst 2D grid of forecast values
     *  @param curve_ref Reference quantiles
     *  @param curve_fcst Forecast quantiles
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @return Calibrated forecasts
    */
    vec2 apply_curve(const vec2& fcst, const vec& curve_ref, const vec& curve_fcst, Extrapolation policy_below, Extrapolation policy_above);

    /** Apply arbitrary calibration curve to 2D forecasts with spatially varying QQ map
     *  @param fcst 2D grid of forecast values
     *  @param curve_ref Reference quantiles (Y, X, Q)
     *  @param curve_fcst Forecast quantiles (Y, X, Q)
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @return Calibrated forecasts
    */
    vec2 apply_curve(const vec2& fcst, const vec3& curve_ref, const vec3& curve_fcst, Extrapolation policy_below, Extrapolation policy_above);

    /** Ensure calibration curve is monotonic, by removing points
     *  @param curve_ref Reference quantiles
     *  @param curve_fcst Forecast quantiles
     *  @param output_fcst New forecast quantiles
     *  @returns New reference quantiles
    */
    vec monotonize_curve(vec curve_ref, vec curve_fcst, vec& output_fcst);

    float get_optimal_threshold(const vec& ref, const vec& fcst, float threshold, Metric metric);

    float calc_score(float a, float b, float c, float d, Metric metric);
    float calc_score(const vec& ref, const vec& fcst, float threshold, Metric metric);
    float calc_score(const vec& ref, const vec& fcst, float threshold, float fthreshold, Metric metric);

    /**@}*/

    /** **************************************
     * @name Downscaling methods
     * Functions that interpolate data from one grid to another
     * ***************************************/ /**@{*/

    /** Generic downscaler for simple methods that don't use information other than the grids
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @param downscaler Downscaling method
      * @return Values on the output grid
    */
    vec2 downscaling(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, Downscaler downscaler);
    vec3 downscaling(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, Downscaler downscaler);
    vec downscaling(const Grid& igrid, const Points& opoints, const vec2& ivalues, Downscaler downscaler);
    vec2 downscaling(const Grid& igrid, const Points& opoints, const vec3& ivalues, Downscaler downscaler);
    // vec downscaling(const Points& ipoints, const Points& opoints, const vec& ivalues, Downscaler downscaler);
    // vec2 downscaling(const Points& ipoints, const Points& opoints, const vec2& ivalues, Downscaler downscaler);
    // vec2 downscaling(const Points& ipoints, const Grid& ogrid, const vec& ivalues, Downscaler downscaler);
    // vec3 downscaling(const Points& ipoints, const Grid& ogrid, const vec2& ivalues, Downscaler downscaler);

    /** Nearest neighbour dowscaling grid to grid
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values on the output grid
    */
    vec2 nearest(const Grid& igrid, const Grid& ogrid, const vec2& ivalues);
    vec3 nearest(const Grid& igrid, const Grid& ogrid, const vec3& ivalues);

    /** Nearest neighbour dowscaling grid to point
      * @param igrid Input grid
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values for the output points
    */
    vec nearest(const Grid& igrid, const Points& opoints, const vec2& ivalues);
    vec2 nearest(const Grid& igrid, const Points& opoints, const vec3& ivalues);

    /** Nearest neighbour dowscaling point to point
      * @param ipoints Input points
      * @param opoints Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values for the output points
    */
    vec nearest(const Points& ipoints, const Points& opoints, const vec& ivalues);
    vec2 nearest(const Points& ipoints, const Points& opoints, const vec2& ivalues);

    /** Nearest neighbour dowscaling point to grid
      * @param ipoints Input points
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values for the output points
    */
    vec2 nearest(const Points& ipoints, const Grid& ogrid, const vec& ivalues);
    vec3 nearest(const Points& ipoints, const Grid& ogrid, const vec2& ivalues);

    /** Nearest neighbour downscaling grid to grid and probability in one
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues 3D vector of values on the input grid (Y, X, E)
      * @param threshold 2D vector of threshold values
      * @param comparison_operator lower than, lower or equal than, greater than, great or equal than
      * @return Values on the output grid
    */
    vec2 downscale_probability(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec2& threshold, const ComparisonOperator& comparison_operator);

    /** Nearest neighbour downscaling grid to grid and threshold and consensus in one
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues_true 3D vector of values on the input grid (Y, X, E)
      * @param ivalues_false 3D vector of values on the input grid (Y, X, E)
      * @param threshold_values 3D vector of values (Y, X, E), which defines the mask array for ivalues after applying the theshold
      * @param threshold 2D vector of threshold values
      * @param comparison_operator lower than, lower or equal than, greater than, great or equal than
      * @param statistic statistic to compute over the ensemble dimension
      * @return Values on the output grid
    */
    vec2 mask_threshold_downscale_consensus(const Grid& igrid, const Grid& ogrid, const vec3& ivalues_true, const vec3& ivalues_false, const vec3& theshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const Statistic& statistic);

    /** Nearest neighbour downscaling grid to grid and threshold and quantile in one
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues_true 3D vector of values on the input grid (Y, X, E)
      * @param ivalues_false 3D vector of values on the input grid (Y, X, E)
      * @param threshold_values 3D vector of values (Y, X, E), which defines the mask array for ivalues after applying the theshold
      * @param threshold 2D vector of threshold values
      * @param comparison_operator lower than, lower or equal than, greater than, great or equal than
      * @param quantile quantile (value between 0-1) to compute
      * @return Values on the output grid
    */
    vec2 mask_threshold_downscale_quantile(const Grid& igrid, const Grid& ogrid, const vec3& ivalues_true, const vec3& ivalues_false, const vec3& theshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const float quantile);

    /** Bilinear downscaling grid to grid
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values on the output grid
    */
    vec2 bilinear(const Grid& igrid, const Grid& ogrid, const vec2& ivalues);
    vec3 bilinear(const Grid& igrid, const Grid& ogrid, const vec3& ivalues);

    /** Bilinear downscaling grid to points
      * @param igrid Input grid
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values for the output points
    */
    vec bilinear(const Grid& igrid, const Points& opoints, const vec2& ivalues);
    vec2 bilinear(const Grid& igrid, const Points& opoints, const vec3& ivalues);

    vec2 simple_gradient(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, float elev_gradient, Downscaler downscaler=Nearest);
    vec3 simple_gradient(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, float elev_gradient, Downscaler downscaler=Nearest);

    vec simple_gradient(const Grid& igrid, const Points& opoints, const vec2& ivalues, float elev_gradient, Downscaler downscaler=Nearest);
    vec2 simple_gradient(const Grid& igrid, const Points& opoints, const vec3& ivalues, float elev_gradient, Downscaler downscaler=Nearest);

    /** Compute Downscale
    *@param igrid input grid
    *@param ogrid output grid
    *@param ivalues values from igrid
    *@param elev_gradient elevation gradient
    *@param laf_gradient land area fraction gradient
    */
    vec2 full_gradient(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient=vec2(), Downscaler downscaler=Nearest);
    vec3 full_gradient(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec3& elev_gradient, const vec3& laf_gradient, Downscaler downscaler=Nearest);
    vec full_gradient(const Grid& igrid, const Points& opoints, const vec2& ivalues, const vec2& elev_gradient, const vec2& laf_gradient, Downscaler downscaler=Nearest);
    vec2 full_gradient(const Grid& igrid, const Points& opoints, const vec3& ivalues, const vec3& elev_gradient, const vec3& laf_gradient, Downscaler downscaler=Nearest);

    /* Elevation and land area fraction downscaling with debug output fields
    */
    vec3 full_gradient_debug(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient=vec2(), Downscaler downscaler=Nearest);

    /** Smart neighbour downscaling grid to grid
      * @param igrid Input grid
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @param num Number of neighbours to average
      * @param structure Structure function for determining similarity
      * @return Values for the output points
    */
    vec2 smart(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, int num, const StructureFunction& structure);
    /**@}*/

    /** **************************************
     * @name Grid calculations
     * Functions that calculate statistics on a grid
     * ***************************************/ /**@{*/

    /** For each point, counts the number of gridpoints within the radius
     *  @param grid Grid
     *  @param points Points
     *  @param radius Radius [m]
     *  @return Number of gridpoints
    */
    vec count(const Grid& grid, const Points& points, float radius);

    /** For each gridpoint, counts the number of gridpoints within the radius
     *  @param igrid Input grid
     *  @param ogrid Output grid
     *  @param radius Radius [m]
     *  @return Number of gridpoints
    */
    vec2 count(const Grid& igrid, const Grid& ogrid, float radius);

    /** For each gridpoint, counts the number of points within the radius
     *  @param grid Grid
     *  @param points Points
     *  @param radius Radius [m]
     *  @return Number of points
    */
    vec2 count(const Points& points, const Grid& grid, float radius);

    /** For each point, counts the number of points within the radius
     *  @param ipoints Input points
     *  @param opoints Output points
     *  @param radius Radius [m]
     *  @return Number of points
    */
    vec count(const Points& ipoints, const Points& opoints, float radius);

    /** For each point, calculates the distance to nearest gridpoint
     *  @param grid Grid
     *  @param points Points
     *  @param num Number of points
     *  @return Distance [m] to nearest gridpoint for each point
    */
    vec distance(const Grid& grid, const Points& points, int num=1);

    /** For each output gridpoint, calculate the distance to nearest input gridpoint
     *  @param grid Grid
     *  @param ogrid Output grid
     *  @param num Number of points
     *  @return Distance [m] to nearest gridpoint for each gridpoint
    */
    vec2 distance(const Grid& igrid, const Grid& ogrid, int num=1);

    /** For each gridpoint, calculates the distance to nearest point
     *  @param points Points
     *  @param grid Grid
     *  @param num Number of points
     *  @return Distance [m] to nearest gridpoint for each point
    */
    vec2 distance(const Points& points, const Grid& grid, int num=1);

    /** For each output point, calculates the distance to nearest input point
     *  @param ipoints Input points
     *  @param opoints Output points
     *  @param num Number of points
     *  @return Distance [m] to nearest point for each point
    */
    vec distance(const Points& ipoints, const Points& opoint, int num=1);

    /** Fill in missing values based on nearby values
      * @param values 2D array of values
      * @return 2D array of values without any missing values
    */
    vec2 fill_missing(const vec2& values);

    /** Aggregate points onto a grid. Writes MV where there are not enough observations.
      * @param grid Grid to aggregate to
      * @param points Points with values
      * @param values Values at points
      * @param radius Circle radius for aggregate points [m]
      * @param min_num Minimum number of points in radius to create a value
      * @param statistic Statistic to compute on points within radius
    */
    vec2 gridding(const Grid& grid, const Points& points, const vec& values, float radius, int min_num, Statistic statistic);

    /** Aggregate points onto a points. Writes MV where there are not enough observations.
      * @param opoints Points to aggregate to
      * @param ipoints Points with values
      * @param values Values at points
      * @param radius Circle radius for aggregate points [m]
      * @param min_num Minimum number of points in radius to create a value
      * @param statistic Statistic to compute on points within radius
    */
    vec gridding(const Points& opoints, const Points& ipoints, const vec& values, float radius, int min_num, Statistic statistic);

    /** Assign each point to nearest neighbour in grid and aggregate values. Writes MV where there are not enough observations.
      * @param grid Grid to aggregate to
      * @param points Points with values
      * @param values Values at points
      * @param min_num Minimum number of points in a gridpoint to create a value
      * @param statistic Statistic to compute on points within gridpoint
    */
    vec2 gridding_nearest(const Grid& grid, const Points& points, const vec& values, int min_num, gridpp::Statistic statistic);

    /** Assign each ipoint to nearest neighbour in opoints and aggregate values. Writes MV where there are not enough observations.
      * @param opoints Points to aggregate to
      * @param ipoints Points with values
      * @param values Values at points
      * @param min_num Minimum number of points in a gridpoint to create a value
      * @param statistic Statistic to compute on points within gridpoint
    */
    vec gridding_nearest(const Points& opoints, const Points& ipoints, const vec& values, int min_num, gridpp::Statistic statistic);

    /** Compute a score for a metric of all points within a radius */
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
     *  @param temperature Temperatures [K]
     *  @param relative_humidity Relative humidities [1]
     *  @returns Dewpoint temperatures [K]
    */
    vec dewpoint(const vec& temperature, const vec& relative_humidity);

    /** Calculate pressure at a new elevation
     *  @param ielev Elevation at start point
     *  @param oelev Elevation at new point
     *  @param ipressure Pressure at start point
     *  @param itemperature Temperature at start point
     *  @return Pressure at new point
     */
    float pressure(float ielev, float oelev, float ipressure, float itemperature=288.15);

    /** Calculate Vector version of pressure calculation
     *  @param ielev Elevations at start point
     *  @param oelev Elevations at new point
     *  @param ipressure Pressures at start point
     *  @param itemperature Temperatures at start point
     *  @return Pressures at new points
     */
    vec pressure(const vec& ielev, const vec& oelev, const vec& ipressure, const vec& itemperature);

    /** Convert Surface Pressure to Sea Level Pressure
     *  @param ps Surface pressure [pa]
     *  @param altitude Station altitude above sea level [m]
     *  @param temperature 2m temperature [K]
     *  @param rh 2m Relative humidity [1]
     *  @param dewpoint 2m Dewpoint Temperature at station [K]
     *  @return Sea Level Pressure [pa]
     */
    float sea_level_pressure(float ps, float altitude, float temperature, float rh=gridpp::MV, float dewpoint=gridpp::MV);

    /** Vector version of Convert Surface Pressure to Sea Level Pressure
     *  @param ps Surface pressures [pa]
     *  @param altitude Station altitudes above sea level [m]
     *  @param temperature 2m temperatures [K]
     *  @param rh 2m Relative humidities [1]
     *  @param dewpoint 2m Dewpoint Temperatures at stations [K]
     *  @return Sea Level Pressure [pa]
     */
    vec sea_level_pressure(const vec& ps, const vec& altitude, const vec& temperature, const vec& rh, const vec& dewpoint);

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

    /** Calculate relative humidity from temperature and dewpoint temperature
     *  @param temperature Temperature [K]
     *  @param dewpoint Dewpoint temperature [K]
     *  @returns Relative humidity [1]
    */
    float relative_humidity(float temperature, float dewpoint);

    /** Vector version of relative humidity calculation
     *  @param temperature Temperatures [K]
     *  @param dewpoint Dewpoint temperatures [K]
     *  @returns Relative humidities [1]
    */
    vec relative_humidity(const vec& temperature, const vec& dewpoint);

    /** Calculate wetbulb temperature from temperature, pressure, and relative humidity
     *  @param temperature Temperature [K]
     *  @param pressure Air pressure [pa]
     *  @param Relative humidity [1]
     *  @returns Wetbulb temperature [K]
    */
    float wetbulb(float temperature, float pressure, float relative_humidity);

    /** Vector version of wetbulb calculation
     *  @param temperature Temperatures [K]
     *  @param pressure Air pressures [pa]
     *  @param Relative humidities [1]
     *  @returns Wetbulb temperatures [K]
    */
    vec wetbulb(const vec& temperature, const vec& pressure, const vec& relative_humidity);

    /** Diagnose wind speed from its components
     *  @param xwind X-component of wind [any unit]
     *  @param ywind Y-component of wind [any unit]
     *  @return Wind speed [any unit]
     * */
    float wind_speed(float xwind, float ywind);

    /** Vector version of wind speed calculation
     *  @param xwind X-components of wind [any unit]
     *  @param ywind Y-components of wind [any unit]
     *  @return Wind speeds [any unit]
     * */
    vec wind_speed(const vec& xwind, const vec& ywind);

    /** Diagnose wind direction from its components. If both xwind and ywind are 0, then direction
     *  is 180
     *  @param xwind X-component of wind [any unit]
     *  @param ywind Y-component of wind [any unit]
     *  @return Wind direction [degrees]
     * */
    float wind_direction(float xwind, float ywind);

    /** Vector version of wind direction calculation
     *  @param xwind X-components of wind [any unit]
     *  @param ywind Y-components of wind [any unit]
     *  @return Wind direction [degrees]
     * */
    vec wind_direction(const vec& xwind, const vec& ywind);

    /**@}*/

    /** ****************************************
     * @name OpenMP settings
     * Functions that configure OpenMP
     * *****************************************/ /**@{*/
    /** Set the number of OpenMP threads to use. Overrides OMP_NUM_THREAD env variable. */
    void set_omp_threads(int num);

    /** Get the number of OpenMP threads currently set */
    int get_omp_threads();

    /** Sets the number of OpenMP threads to 1 if OMP_NUM_THREADS undefined */
    void initialize_omp();

    /** ****************************************
     * @name Utilities
     * Helper functions
     * *****************************************/ /**@{*/
    // ivec regression(const vec& x, const vec& y);

    /** Set the verbosity of debug messages. Use 0 for no messages. */
    void set_debug_level(int level);

    static int _debug_level = 0;

    /** Get the currently set level of debug messages */
    int get_debug_level();

    /** Convert name of a statistic enum */
    Statistic get_statistic(std::string name);

    /** The gridpp version
     * @return The gridpp version
    */
    std::string version();

    double clock();
    void debug(std::string string);
    void warning(std::string string);
    void error(std::string string);
    void future_deprecation_warning(std::string function, std::string other="");
    bool is_valid(float value);
    float calc_statistic(const vec& array, Statistic statistic);
    float calc_quantile(const vec& array, float quantile);
    vec calc_statistic(const vec2& array, Statistic statistic);
    vec calc_quantile(const vec2& array, float quantile=MV);

    /** Compute quantile with 2D varying quantile
      * @param array Input array with dimensions (T, Y, X)
      * @param quantile Quantile array with dimensions (Y, X)
      * @return Extracted quantile with dimensions (Y, X)
    */
    vec2 calc_quantile(const vec3& array, const vec2& quantile);
    int num_missing_values(const vec2& iArray);

    /** Find the index in a vector that is equal or just below a value
     *  @param iX Lookup value
     *  @param iValues Lookup vector. Must be sorted.
     *  @return The index into iValues that is equal or just below iX
    */
    int get_lower_index(float iX, const vec& iValues);

    /** Find the index in a vector that is equal or just above a value
     *  @param iX Lookup value
     *  @param iValues Lookup vector. Must be sorted.
     *  @return The index into iValues that is equal or just above iX
    */
    int get_upper_index(float iX, const vec& iValues);

    /** Piecewise linear interpolation
     *  If x is outside the range of iX, then the min/max value of iY is used. If there are multiple
     *  identical x values, then the average of the y values at each end of the x-interval that is
     *  constant is used. The exception is if the constant interval is on one (only) of the edges of the
     *  interpolation curve. In that case, the y-value at the end of the interval further away from
     *  the boundary of the curve is used.
     *  @param x Interpolation to this point
     *  @param iX Vector of x-axis values. Vector must be sorted.
     *  @param iY Vector of y-axis values corresponding to iX.
     *  @return Y value corresponding to x
    */
    float interpolate(float x, const vec& iX, const vec& iY);

    /** Piecewise linear interpolation
     *  Same as above, but for a vector of points to interpolate for
     *  @param x Interpolation to these points
     *  @param iX Vector of x-axis values. Vector must be sorted.
     *  @param iY Vector of y-axis values corresponding to iX.
     *  @return Y values corresponding to x
    */
    vec interpolate(const vec& x, const vec& iX, const vec& iY);

    /** Initialize a vector of size Y, X, with a given value */
    ivec2 init_ivec2(int Y, int X, int value);
    vec2 init_vec2(int Y, int X, float value=MV);

    /** Initialize a vector of size Y, X, E, with a given value */
    ivec3 init_ivec3(int Y, int X, int E, int value);
    vec3 init_vec3(int Y, int X, int E, float value=MV);

    /** Get reasonably spaced quantiles from a vector of values, ignoring duplicate values
      *  but including the first number after duplicated values. Include the lowest and highest
      *  values.
      *  @param values vector of values (unsorted, and no invalid values)
      *  @param num number of thresholds to get
    */
    vec calc_even_quantiles(const vec& values, int num);

    /** Computes gradients based on values in neighbourhood
     *  @param grid Grid
     *  @param base Dependent variable. Missing values are not used.
     *  @param values Independent variable. Missing values are not used.
     *  @param radius Neighbourhood radius in number of gridpoints
     *  @param min_nim Minimum number of points required to compute gradient
     *  @param min_range Minimum range of base to compute gradient
     *  @param default_gradient Use this gradient if minimum number is not met
    */
    vec2 calc_gradient(const vec2& base, const vec2& values, GradientType gradient_type, int halfwidth, int min_num=2, float min_range=gridpp::MV, float default_gradient=0);
    /** Find suitable value in neighbourhood based on a search criteria. If search value is within a
    * criteria range, then the most suitable point is used. This is the nearest value of any point
    * within the search_target range; or if no point fulfills this, the point with the highest
    * search value.
    * @param base base values (e.g elevation)
    * @param values values to compute gradients for (e.g. temperatures)
    * @param gradient_type what gradient type to compute
    * @param halfwidth neighbourhood halfwidth to compute gradient for
    * @param min_num minimum number of valid points needed to compute gradient
    * @param min_range minimum range of base values to compute gradient (for LinearRegression, this is the standard deviation of values
    * @param default_gradient The gradient to use when a gradient cannot be computed
    */
    vec2 neighbourhood_search(const vec2& array, const vec2& search_array, int halfwidth, float search_target_min, float search_target_max, float search_delta, const ivec2& apply_array=ivec2());

    /** Compute window statistics
    *  @param array input array with dimensions (case, time)
    *  @param length window length in number of timesteps
    *  @param statistic statistic to apply to window
    *  @param before if true, make the window end at the particular time. If false, centre it.
    *  @param keep_missing if true, window value will be missing if one or more values in window are missing
    *  @param missing_edges if true put missing values at the edges, where window overshoots the edge
    */

    vec2 window(const vec2& array, int length, gridpp::Statistic statistic, bool before=false, bool keep_missing=false, bool missing_edges=true);

    /** Check if the grid is the same size as the 2D vector. If True, they are compatible, if false
    * they are incompatible */
    bool compatible_size(const Grid& grid, const vec2& v);
    bool compatible_size(const Grid& grid, const vec3& v);
    bool compatible_size(const Points& points, const vec& v);
    bool compatible_size(const Points& points, const vec2& v);
    bool compatible_size(const vec2& a, const vec2& b);
    bool compatible_size(const vec2& a, const vec3& b);
    bool compatible_size(const vec3& a, const vec3& b);

    /** Checks if a point is located inside a rectangle formed by 4 points. The 4 points must be
      * provided in an order that draws out a rectangle (either clockwise or counter-clockwise)
      * @param A: A point in the rectangle
      * @param B: A point in the rectangle
      * @param C: A point in the rectangle
      * @param D: A point in the rectangle
      * @param m: The point to test if it is inside
      * @return True if the point is inside, False otherwise
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

    /** Default value used to fill array in SWIG testing functions. Not useful for any other purpose. */
    static const float swig_default_value = -1;

    /**@}*/

    /** Represents a single point in some coordinate system */
    class Point {
        public:
            /** Constructor
              * @param lat: Latitude coordinate
              * @param lon: Longitude coordinate
              * @param elev: Elevation
              * @param laf: Land area fraction (between 0 and 1)
              * @param type: Coordinate type for lat and lon
            */
            Point(float lat, float lon, float elev=MV, float laf=MV, CoordinateType type=Geodetic);
            float lat;
            float lon;
            float elev;
            float laf;
            CoordinateType type;
    };

    /** Helper class for Grid and Points */
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


            /** Convert lat/lons to 3D cartesian coordinates with the centre of the earth as the origin
             *  @param lats vector of latitudes [deg]
             *  @param lons vector of longitudes [deg]
             *  @param x_coords vector of x-coordinates [m]
             *  @param y_coords vector of y-coordinates [m]
             *  @param z_coords vector of z-coordinates [m]
             * */
            bool convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords) const;

            /** Same as above, but convert a single lat/lon to 3D cartesian coordinates
             *  @param lat latitude [deg]
             *  @param lon longitude [deg]
             *  @param x_coord x-coordinate [m]
             *  @param y_coord y-coordinate [m]
             *  @param z_coord z-coordinate [m]
             * */
            bool convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord) const;
            static float deg2rad(float deg);
            static float rad2deg(float deg);
            static float calc_distance(float lat1, float lon1, float lat2, float lon2, CoordinateType type=Geodetic);
            static float calc_distance(float x0, float y0, float z0, float x1, float y1, float z1);
            static float calc_distance_fast(float lat1, float lon1, float lat2, float lon2, CoordinateType type=Geodetic);
            static float calc_distance_fast(const Point& p1, const Point& p2);
            vec get_lats() const;
            vec get_lons() const;
            int size() const;
            CoordinateType get_coordinate_type() const;
        protected:
            typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
            typedef std::pair<point, unsigned> value;
            typedef boost::geometry::model::box<point> box;
            boost::geometry::index::rtree< value, boost::geometry::index::quadratic<16> > mTree;
            vec mLats;
            vec mLons;
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
            // Checks that a lat-coordinate is valid (based on the coordinate type)
            bool check_lat(float lat) const;
            // Checks that a lon-coordinate is valid (based on the coordinate type)
            bool check_lon(float lon) const;
    };

    /** Represents a vector of locations and their metadata */
    class Points  {
        public:
            Points();
            /** Initialize a new grid
             *  @param lats: vector of latitudes [degrees]
             *  @param lons: vector of longitudes [degrees]
             *  @param elevs: vector of elevations [m]
             *  @param lafs: vector of land area fractions [1]
             *  @param type: Coordinate type
            */
            Points(vec lats, vec lons, vec elevs=vec(), vec lafs=vec(), CoordinateType type=Geodetic);
            Points(KDTree tree, vec elevs=vec(), vec lafs=vec());
            Points& operator=(Points other);
            Points(const Points& other);
            // Returns -1 if there are no neighbours
            int get_nearest_neighbour(float lat, float lon, bool include_match=true) const;
            ivec get_neighbours(float lat, float lon, float radius, bool include_match=true) const;
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match=true) const;
            int get_num_neighbours(float lat, float lon, float radius, bool include_match=true) const;
            ivec get_closest_neighbours(float lat, float lon, int num, bool include_match=true) const;

            vec get_lats() const;
            vec get_lons() const;
            vec get_elevs() const;
            vec get_lafs() const;
            int size() const;
            ivec get_in_domain_indices(const Grid& grid) const;
            Points get_in_domain(const Grid& grid) const;
            CoordinateType get_coordinate_type() const;
            Point get_point(int index) const;
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
             *  @param lats: 2D vector of latitudes [degrees]
             *  @param lons: 2D vector of longitudes [degrees]
             *  @param elevs: 2D vector of elevations [m]
             *  @param lafs: 2D vector of land area fractions [1]
             *  @param type: Coordinate type
            */
            Grid(vec2 lats, vec2 lons, vec2 elevs=vec2(), vec2 lafs=vec2(), CoordinateType type=Geodetic);
            ivec get_nearest_neighbour(float lat, float lon, bool include_match=true) const;
            ivec2 get_neighbours(float lat, float lon, float radius, bool include_match=true) const;
            ivec2 get_neighbours_with_distance(float lat, float lon, float radius, vec& distances, bool include_match=true) const;
            int get_num_neighbours(float lat, float lon, float radius, bool include_match=true) const;
            ivec2 get_closest_neighbours(float lat, float lon, int num, bool include_match=true) const;

            bool get_box(float lat, float lon, int& Y1_out, int& X1_out, int& Y2_out, int& X2_out) const;

            /** Convert grid to a vector of points */
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
            not_implemented_exception() : std::logic_error("Function not yet implemented") { };
    };

    /** Covariance structure function */
    class StructureFunction {
        public:
            StructureFunction(float localization_distance=0);
            /** Correlation between two points */
            virtual float corr(const Point& p1, const Point& p2) const = 0;
            virtual float corr_background(const Point& p1, const Point& p2) const;
            /** Maximum distance for which an observation can have an impact (localization)
              * @return Distance [m]
            */
            virtual float localization_distance(const Point& p) const;
            virtual StructureFunction* clone() const = 0;
            static const float default_min_rho;
        protected:
            /** Barnes correlation function
              * @param dist Distance between points. Same units as 'length'
              * @param length Length scale
            */
            float barnes_rho(float dist, float length) const;

            /** Cressman correlation function
              * @param dist Distance between points. Same units as 'length'
              * @param length Length scale
            */
            float cressman_rho(float dist, float length) const;
            float m_localization_distance;
    };
    class MultipleStructure: public StructureFunction {
        public:
            /** Exponential structure function
              * @param h: Horizontal decorrelation length >=0 [m]
              * @param v: Vertical decorrelation length >=0 [m]. If 0, disable decorrelation.
              * @param w: Land/sea decorrelation length >=0 [1]. If 0, disable decorrelation.
              * @param hmax: Truncate horizontal correlation beyond this length [m]. If undefined, 3.64 * h.
            */
            MultipleStructure(const StructureFunction& structure_h, const StructureFunction& structure_v, const StructureFunction& structure_w);
            float corr(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
            float localization_distance(const Point& p) const;
        private:
            StructureFunction* m_structure_h;
            StructureFunction* m_structure_v;
            StructureFunction* m_structure_w;
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
              * @param h: Horizontal decorrelation length field >=0, same size as grid [m]
              * @param v: Vertical decorrelation length field >=0 [m]. Set all to 0 to disable decorrelation.
              * @param w: Land/sea decorrelation length field >=0 [1]. Set all to 0 to disable decorrelation.
              * @param min_rho: Truncate horizontal correlation when rho less than this value [m].
            */
            BarnesStructure(Grid grid, vec2 h, vec2 v, vec2 w, float min_rho=StructureFunction::default_min_rho);
            float corr(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
            float localization_distance(const Point& p) const;
        private:
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
            StructureFunction* clone() const;
        private:
            float mH;
            float mV;
            float mW;
    };

    class CrossValidation: public StructureFunction {
        public:
            /** Structure function for performing cross validation experiments
              * @param dist: Force background-to-obs correlation to 0 for points within
              *   this distance [m]. If MV, disable this.
            */
            CrossValidation(StructureFunction& structure, float dist=MV);
            float corr(const Point& p1, const Point& p2) const;
            float corr_background(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
            float localization_distance(const Point& p) const;
        private:
            StructureFunction* m_structure;
            float m_dist;
    };

    class Transform {
        public:
            // Note these cannot be pure virtual, otherwise SWIG does not expose
            // the vector functions (with the same name) in python. Therefore, make sure these
            // functions are overloaded in the subclass implementation
            virtual float forward(float value) const;
            virtual float backward(float value) const;

            vec forward(const vec& input) const;
            vec backward(const vec& input) const;
            vec2 forward(const vec2& input) const;
            vec2 backward(const vec2& input) const;
            vec3 forward(const vec3& input) const;
            vec3 backward(const vec3& input) const;
    };
    class Identity : public Transform {
        public:
            // SWIG requires these "using" statements to enable the vectorized versions in the
            // subclasses
            using Transform::forward;
            using Transform::backward;
            float forward(float value) const;
            float backward(float value) const;
    };
    class Log : public Transform {
        public:
            using Transform::forward;
            using Transform::backward;
            float forward(float value) const;
            float backward(float value) const;
    };
    class BoxCox : public Transform {
        public:
            BoxCox(float threshold);
            using Transform::forward;
            using Transform::backward;
            float forward(float value) const;
            float backward(float value) const;
        private:
            float mThreshold;
    };
    class Gamma : public Transform {
        /** Transforms values to cdf from a gamma distribution and subsequantly extracts the
         *  cdf from a standard normal distribution.
        */
        public:
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
