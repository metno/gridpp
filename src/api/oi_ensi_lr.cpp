#include "gridpp.h"
#include <math.h>
#include <algorithm>
#include <armadillo>
#include <assert.h>
#include <exception>

using namespace gridpp;

namespace {
    typedef arma::mat mattype;
    typedef arma::vec vectype;
    typedef arma::cx_mat cxtype;

    void check_vec(vec2 input, int Y, int X);
    void check_vec(vec input, int S);

    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };
}

template<class Matrix>
void print_matrix(Matrix matrix) {
       matrix.print(std::cout);
}

template void print_matrix< ::mattype>(::mattype matrix);
template void print_matrix< ::cxtype>(::cxtype matrix);

/*vec3 gridpp::optimal_interpolation_ensi_lr(const gridpp::Grid& bgrid,
        const vec3& background_l, 
        const vec3& background_L, 
        const gridpp::Points& points,
        const vec2& pobs,
        const vec2& pbackground_r, 
        const vec2& pbackground_R,
        const gridpp::StructureFunction& structure,
        float var_ratios_or, 
        float std_ratios_lr,
        float weight, 
        int max_points,
        bool allow_extrapolation) {
    double s_time = gridpp::clock();

    // Check input data
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");

    int nS = points.size();
    if(nS == 0)
        return background_l;

    int nY = bgrid.size()[0];
    int nX = bgrid.size()[1];

    if(nY == 0 || nX == 0) {
        std::stringstream ss;
        ss << "Grid size (" << nY << "," << nX << ") cannot be zero";
        throw std::invalid_argument(ss.str());
    }

    if(bgrid.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background grid and observations points must be of same coordinate type (lat/lon or x/y)");
    }
    // Check ensembles have consistent sizes
    int nE = background_l[0][0].size();
    if(background_l.size() != nY || background_l[0].size() != nX) {
        std::stringstream ss;
        ss << "Input left field (" << background_l.size() << "," << background_l[0].size() << "," << background_l[0][0].size() << ") is not the same size as the grid (" << nY << "," << nX << "," << nE << ")";
        throw std::invalid_argument(ss.str());
    }
    if(background_L.size() != nY || background_L[0].size() != nX || background_L[0][0].size() != nE) {
        std::stringstream ss;
        ss << "Input LEFT field (" << background_L.size() << "," << background_L[0].size() << "," << background_L[0][0].size() << ") is not the same size as the grid (" << nY << "," << nX << "," << nE << ")";
        throw std::invalid_argument(ss.str());
    }
    if(pbackground_r.size() != nS || pbackground_r[0].size() != nE) {
        std::stringstream ss;
        ss << "Input right field at observation location (" << pbackground_r.size() << "," << pbackground_r[0].size() << ") and points (" << nS << "," << nE << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }

    if(pbackground_R.size() != nS || pbackground_R[0].size() != nE) {
        std::stringstream ss;
        ss << "Input RIGHT field at observation location (" << pbackground_R.size() << "," << pbackground_R[0].size() << ") and points (" << nS << "," << nE << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    // Check observations have consistent size
    if(pobs.size() != nS || pobs[0].size() != nE) {
        std::stringstream ss;
        ss << "Observations (" << pobs.size() << "," << pobs[0].size() << ") and points (" << nS << "," << nE << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    
    gridpp::Points bpoints = bgrid.to_points();
    vec2 background_l1 = gridpp::init_vec2(nY * nX, nE);
    vec2 background_L1 = gridpp::init_vec2(nY * nX, nE);
    int count = 0;
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            for(int e = 0; e < nE; e++) {
                background_l1[count][e] = background_l[y][x][e];
                background_L1[count][e] = background_L[y][x][e];
            }
            count++;
        }
    }
    vec2 output1 = optimal_interpolation_ensi_lr(bpoints, background_l1, background_L1, points, pobs, pbackground_r, pbackground_R, structure, var_ratios_or, std_ratios_lr, weight, max_points, allow_extrapolation);
    vec3 output = gridpp::init_vec3(nY, nX, nE);
    count = 0;
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            for(int e = 0; e < nE; e++) {
                output[y][x][e] = output1[count][e];
            }
            count++;
        }
    }
    return output;
} */
/* command to fix the gridpp.R file:
for i in {1..50}; do sed -i "s/all(sapply(argv\[\[$i\]\] , is.integer) || sapply(argv\[\[$i\]\], is.numeric))/all(sapply(argv\[\[$i\]\] , is.integer) | sapply(argv\[\[$i\]\], is.numeric))/g" gridpp.R; done */
vec2 gridpp::R_optimal_interpolation_ensi_lr(const gridpp::Points& bpoints,
        const vec2& background_l, 
        const vec2& background_L, 
        const gridpp::Points& points,
        const vec2& pobs,  
        const vec2& pbackground_r, 
        const vec2& pbackground_R,
/*        const gridpp::StructureFunction& structure, it creates problem with R bindings */
        int which_structfun,
        float dh,
        float dz,
        float dw,
        float var_ratios_or,
        float std_ratios_lr,
        float weight,
        int max_points,
        bool allow_extrapolation) {
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");
    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background and observations points must be of same coorindate type (lat/lon or x/y)");
    }
    if(background_l.size() != bpoints.size())
        throw std::invalid_argument("Input left field is not the same size as the grid");
    if(background_L.size() != bpoints.size())
        throw std::invalid_argument("Input LEFT field is not the same size as the grid");
    if(pobs.size() != points.size())
        throw std::invalid_argument("Observations and points exception mismatch");
    if(pbackground_r.size() != points.size())
        throw std::invalid_argument("Background rigth and points size mismatch");
    if(pbackground_R.size() != points.size())
        throw std::invalid_argument("Background RIGTH and points size mismatch");

    float hmax = 7 * dh;
    float default_min_std = 0.0013;
/*
    BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax)t;
    if(which_structfun == 0) {
        BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax);
    }
    else if(which_structfun == 1) {
        MixAStructure structure = MixAStructure( dh, dz, dw, hmax);
    } 
    else {
        BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax);
    } */

    int nS = points.size();
    if(nS == 0)
        return background_l;

    int mY = -1;  // Write debug information for this station index
    bool diagnose = false;

    int nY = background_l.size();
    int nEns = background_l[0].size();

    // Prepare output matrix
    float missing_value = -99999.999;
/*    vec2 output = gridpp::init_vec2(nY, nEns, missing_value); */
    vec2 output = background_l;

    vec blats = bpoints.get_lats();
    vec blons = bpoints.get_lons();
    vec belevs = bpoints.get_elevs();
    vec blafs = bpoints.get_lafs();

    vec plats = points.get_lats();
    vec plons = points.get_lons();
    vec pelevs = points.get_elevs();
    vec plafs = points.get_lafs();


    // Create all objects of type Point (to save time on creation later)
    std::vector<Point> point_vec;
    point_vec.reserve(nS);
    for(int i = 0; i < nS; i++) {
        point_vec.push_back(points.get_point(i));
    }

    // Calculate number of valid members
    // NOTE: all backgrounds must be simultaneously valid
    int nValidEns = 0;
    ivec validEns;
    for(int e = 0; e < nEns; e++) {
        int numInvalid = 0;
        for(int y = 0; y < nY; y++) {
            float value_l = background_l[y][e];
            float value_L = background_L[y][e];
            if(!gridpp::is_valid(value_l) || !gridpp::is_valid(value_L))
                numInvalid++;
        }
        for(int i = 0; i < nS; i++) {
            float value_r = pbackground_r[i][e];
            float value_R = pbackground_R[i][e];
            if(!gridpp::is_valid(value_r) || !gridpp::is_valid(value_R))
                numInvalid++;
        }
        if(numInvalid == 0) {
            validEns.push_back(e);
            nValidEns++;
        }
    }
    if(nValidEns == 0)
        return background_l;

    // gZ_R(nY, nValidEns): used to compute ensemble-based background correlations i) between yth gridpoint and observations ii) among observations
    vec2 gZ_R = gridpp::init_vec2(nY, nValidEns); // useful to compute dynamical correlations
    for(int i = 0; i < nS; i++) {
        vec pbackgroundValid_R(nValidEns);
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            pbackgroundValid_R[e] = pbackground_R[i][ei];
        }
        float mean = gridpp::calc_statistic(pbackgroundValid_R, gridpp::Mean);
        float std = gridpp::calc_statistic(pbackgroundValid_R, gridpp::Std);
        if(!gridpp::is_valid(mean) || !gridpp::is_valid(std)) {
            for(int e = 0; e < nValidEns; e++) 
                gZ_R[i][e] = 0;
        }
        else {
            if(std <= default_min_std) {
                for(int e = 0; e < nValidEns; e++) 
                    gZ_R[i][e] = 0;
            }
            else {
                for(int e = 0; e < nValidEns; e++) 
                    gZ_R[i][e] = 1 / sqrt(nValidEns-1) * (pbackgroundValid_R[e] - mean) / std;
            }
        }
    } // end loop over all observations

    // This causes segmentation fault when building the gridpp pypi package
    // 1) Tested removing num_condition_warning and num_real_part_warning from the parallel loop
    //    but it doesnt seem to help
    // #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float lat = blats[y];
        float lon = blons[y];
        float elev = belevs[y];
        float laf = blafs[y];
        Point p1 = bpoints.get_point(y);
        /* float localizationRadius = structure.localization_distance(p1); */
        float localizationRadius = 0;
        if(which_structfun == 0) {
/*            BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax);*/
            BarnesStructure structure = BarnesStructure( dh, dz, dw);
            localizationRadius = structure.localization_distance(p1);
        }
        else if(which_structfun == 1) {
/*            MixAStructure structure = MixAStructure( dh, dz, dw, hmax); */
            MixAStructure structure = MixAStructure( dh, dz, dw);
            localizationRadius = structure.localization_distance(p1);
        }
 
        // Create list of locations for this gridpoint
        ivec lLocIndices0 = points.get_neighbours(lat, lon, localizationRadius);
        if(lLocIndices0.size() == 0) {
            // If we have too few observations though, then use the background
            continue;
        }
        ivec lLocIndices;
        lLocIndices.reserve(lLocIndices0.size());
        std::vector<std::pair<float,int> > lRhos0;

        // Calculate gridpoint to observation rhos
        lRhos0.reserve(lLocIndices0.size());
        std::vector<Point> p2;
        p2.reserve(lLocIndices0.size());
        for(int i = 0; i < lLocIndices0.size(); i++) {
            int index = lLocIndices0[i];
            p2.push_back(point_vec[index]);
        }
        vec rhos(lLocIndices0.size());
        if(which_structfun == 0) {
/*            BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax); */
            BarnesStructure structure = BarnesStructure( dh, dz, dw);
            rhos = structure.corr_background(p1, p2);
        }
        else if(which_structfun == 1) {
/*            MixAStructure structure = MixAStructure( dh, dz, dw, hmax); */
            MixAStructure structure = MixAStructure( dh, dz, dw);
            rhos = structure.corr_background(p1, p2);
        } 
/*        vec rhos = structure.corr_background(p1, p2); */
        for(int i = 0; i < lLocIndices0.size(); i++) {
            int index = lLocIndices0[i];
            if(gridpp::is_valid(pobs[index][0])) {
                if(rhos[i] > 0) {
                    lRhos0.push_back(std::pair<float,int>(rhos[i], i));
                }
            }
        }

        // Make sure we don't use too many observations
        arma::vec lRhos;
        if(max_points > 0 && lRhos0.size() > max_points) {
            // If sorting is enabled and we have too many locations, then only keep the best ones based on rho.
            // Otherwise, just use the last locations added
            lRhos = arma::vec(max_points);
            std::sort(lRhos0.begin(), lRhos0.end(), ::sort_pair_first<float,int>());
            for(int i = 0; i < max_points; i++) {
                // The best values start at the end of the array
                int index = lRhos0[lRhos0.size() - 1 - i].second;
                lLocIndices.push_back(lLocIndices0[index]);
                lRhos(i) = lRhos0[lRhos0.size() - 1 - i].first;
            }
        }
        else {
            lRhos = arma::vec(lRhos0.size());
            for(int i = 0; i < lRhos0.size(); i++) {
                int index = lRhos0[i].second;
                lLocIndices.push_back(lLocIndices0[index]);
                lRhos(i) = lRhos0[i].first;
            }
        }

        int lS = lLocIndices.size();
        if(lS == 0) {
            // If we have too few observations though, then use the background
            continue;
        }

        vectype lElevs(lS);
        vectype lLafs(lS);
        for(int i = 0; i < lLocIndices.size(); i++) {
            int index = lLocIndices[i];
            lElevs[i] = pelevs[index];
            lLafs[i] = plafs[index];
        }

        // lX_L: used to compute ensemble-based background correlations between yth gridpoint and observations
        mattype lX_L(1, nValidEns, arma::fill::zeros);
        vec backgroundValid_L(nValidEns);
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            backgroundValid_L[e] = background_L[y][ei];
        }
        float mean = gridpp::calc_statistic(backgroundValid_L, gridpp::Mean);
        float std = gridpp::calc_statistic(backgroundValid_L, gridpp::Std);
        if(gridpp::is_valid(mean) && gridpp::is_valid(std) && std > default_min_std) {
            for(int e = 0; e < nValidEns; e++) 
                lX_L(0,e) = 1 / sqrt(nValidEns-1) * (backgroundValid_L[e] - mean) / std;
        }
        // lZ_R: used to compute ensemble-based background correlations i) between yth gridpoint and observations ii) among observations
        mattype lZ_R(lS, nValidEns);
        // Innovation: Observation - Background_r
        mattype lInnov(lS, nValidEns, arma::fill::zeros);
        // lR_dd: Observation error correlation matrix
        mattype lR_dd(lS, lS, arma::fill::zeros);
        // lLoc1D: localization for ensemble-based background correlations between yth gridpoint and observations
        mattype lLoc1D(1, lS, arma::fill::zeros); 
        // lLoc2D: localization for ensemble-based background correlations among observations
        mattype lLoc2D(lS, lS, arma::fill::zeros); 
        for(int i = 0; i < lS; i++) {
            // lR_dd: diagonal observation error correlation matrix
            lR_dd(i, i) = var_ratios_or;
            // lLoc1D between yth gridpoint and ith observation computed before
            lLoc1D(0, i) = lRhos(i);
            // compute lZ_R and lInnov
            int index = lLocIndices[i];
            for(int e = 0; e < nValidEns; e++) {
                lZ_R(i, e) = gZ_R[index][e];
                int ei = validEns[e];
                lInnov(i, ei) = pobs[index][ei] - pbackground_r[index][ei];
            }
            // compute lLoc2D
            Point p1 = point_vec[index];
            std::vector<Point> p2(lS, Point(0, 0));
            for(int j = 0; j < lS; j++) {
                int index_j = lLocIndices[j];
                p2[j] = point_vec[index_j];
            }
            vec corr(lS);
            if(which_structfun == 0) {
/*                BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax); */
                BarnesStructure structure = BarnesStructure( dh, dz, dw);
                corr = structure.corr_background(p1, p2);
            }
            else if(which_structfun == 1) {
/*                MixAStructure structure = MixAStructure( dh, dz, dw, hmax); */
                MixAStructure structure = MixAStructure( dh, dz, dw);
                corr = structure.corr_background(p1, p2);
            } 
            /* vec corr = structure.corr(p1, p2); */
            for(int j = 0; j < lS; j++) 
                lLoc2D(i, j) = corr[j];
        } // end loop over closer observations
        // lr_lr(1, lS): ensemble-based background correlations between yth gridpoint and observations
        //  note: correlation between the LEFT background and the RIGHT background (..._lr)
        mattype lr_lr = lLoc1D % (lX_L * lZ_R.t());
        // lR_rr(lS, lS): ensemble-based background correlations among observations
        //  note: correlation between the RIGTH background and the RIGHT background (..._rr)
        mattype lR_rr = lLoc2D % (lZ_R * lZ_R.t());
        // lK(1, lS): Kalman gain
        mattype lK = lr_lr * arma::inv(lR_rr + lR_dd);
        // dx(1, nValidEns): analysis increment 
        mattype dx = std_ratios_lr * weight * (lK * lInnov);
        ///////////////////////////////
        // Anti-extrapolation filter //
        ///////////////////////////////
        if(!allow_extrapolation) {
            arma::rowvec maxIncEns = arma::max(lInnov);
            arma::rowvec minIncEns = arma::min(lInnov);

            for(int e = 0; e < nValidEns; e++) {
                float increment = dx[e];
                float maxInc = maxIncEns[e];
                float minInc = minIncEns[e];
                if(maxInc > 0 && increment > maxInc) {
                   increment = maxInc;
                }
                else if(maxInc < 0 && increment > 0) {
/*                   increment = maxInc; */
                   increment = 0;
                }
                else if(minInc < 0 && increment < minInc) {
                   increment = minInc;
                }
                else if(minInc > 0 && increment < 0) {
/*                   increment = minInc; */
                   increment = 0;
                }
                dx[e] = increment;
            }
        }

        // Compute the analysis as the updated background_l
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            output[y][ei] = background_l[y][ei] + dx[e];
        }
        // debug
/*        for(int i = 0; i < lS; i++) {
            // compute lZ_R and lInnov
            int index = lLocIndices[i];
            std::cout << i << " backg_r obs " << pbackground_r[index][0] << " " << pobs[index][0] << std::endl;
        } // end loop over closer observations
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            std::cout << ei << " backg_l analysis " << background_l[y][ei] << " " << output[y][ei] << std::endl;
        } */

    } // end loop over gridpoint 
    return output;
} // end of R_optimal_interpolation_ensi_lr

/* command to fix the gridpp.R file:
for i in {1..50}; do sed -i "s/all(sapply(argv\[\[$i\]\] , is.integer) || sapply(argv\[\[$i\]\], is.numeric))/all(sapply(argv\[\[$i\]\] , is.integer) | sapply(argv\[\[$i\]\], is.numeric))/g" gridpp.R; done */
vec2 gridpp::R_optimal_interpolation_ensi_staticcorr_lr(const gridpp::Points& bpoints,
        const vec2& background_l, 
        const gridpp::Points& points,
        const vec2& pobs,  
        const vec2& pbackground_r, 
/*        const gridpp::StructureFunction& structure, it creates problem with R bindings */
        int which_structfun,
        float dh,
        float dz,
        float dw,
        float var_ratios_or,
        float std_ratios_lr,
        float weight,
        int max_points,
        bool allow_extrapolation) {
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");
    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background and observations points must be of same coorindate type (lat/lon or x/y)");
    }
    if(background_l.size() != bpoints.size())
        throw std::invalid_argument("Input left field is not the same size as the grid");
    if(pobs.size() != points.size())
        throw std::invalid_argument("Observations and points exception mismatch");
    if(pbackground_r.size() != points.size())
        throw std::invalid_argument("Background rigth and points size mismatch");

    float hmax = 7 * dh;
    float default_min_std = 0.0013;

    int nS = points.size();
    if(nS == 0)
        return background_l;

    int mY = -1;  // Write debug information for this station index
    bool diagnose = false;

    int nY = background_l.size();
    int nEns = background_l[0].size();

    // Prepare output matrix
    float missing_value = -99999.999;
    vec2 output = background_l;

    vec blats = bpoints.get_lats();
    vec blons = bpoints.get_lons();
    vec belevs = bpoints.get_elevs();
    vec blafs = bpoints.get_lafs();

    vec plats = points.get_lats();
    vec plons = points.get_lons();
    vec pelevs = points.get_elevs();
    vec plafs = points.get_lafs();

    // Create all objects of type Point (to save time on creation later)
    std::vector<Point> point_vec;
    point_vec.reserve(nS);
    for(int i = 0; i < nS; i++) {
        point_vec.push_back(points.get_point(i));
    }

    // Calculate number of valid members
    // NOTE: all backgrounds must be simultaneously valid
    int nValidEns = 0;
    ivec validEns;
    for(int e = 0; e < nEns; e++) {
        int numInvalid = 0;
        for(int y = 0; y < nY; y++) {
            float value_l = background_l[y][e];
            if(!gridpp::is_valid(value_l))
                numInvalid++;
        }
        for(int i = 0; i < nS; i++) {
            float value_r = pbackground_r[i][e];
            if(!gridpp::is_valid(value_r))
                numInvalid++;
        }
        if(numInvalid == 0) {
            validEns.push_back(e);
            nValidEns++;
        }
    }
    if(nValidEns == 0)
        return background_l;

    // This causes segmentation fault when building the gridpp pypi package
    // 1) Tested removing num_condition_warning and num_real_part_warning from the parallel loop
    //    but it doesnt seem to help
    // #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float lat = blats[y];
        float lon = blons[y];
        float elev = belevs[y];
        float laf = blafs[y];
        Point p1 = bpoints.get_point(y);
        /* float localizationRadius = structure.localization_distance(p1); */
        float localizationRadius = 0;
        if(which_structfun == 0) {
/*            BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax); */
            BarnesStructure structure = BarnesStructure( dh, dz, dw);
            localizationRadius = structure.localization_distance(p1);
        }
        else if(which_structfun == 1) {
/*            MixAStructure structure = MixAStructure( dh, dz, dw, hmax); */
            MixAStructure structure = MixAStructure( dh, dz, dw);
            localizationRadius = structure.localization_distance(p1);
        }
 
        // Create list of locations for this gridpoint
        ivec lLocIndices0 = points.get_neighbours(lat, lon, localizationRadius);
        if(lLocIndices0.size() == 0) {
            // If we have too few observations though, then use the background
            continue;
        }
        ivec lLocIndices;
        lLocIndices.reserve(lLocIndices0.size());
        std::vector<std::pair<float,int> > lRhos0;

        // Calculate gridpoint to observation rhos
        lRhos0.reserve(lLocIndices0.size());
        std::vector<Point> p2;
        p2.reserve(lLocIndices0.size());
        for(int i = 0; i < lLocIndices0.size(); i++) {
            int index = lLocIndices0[i];
            p2.push_back(point_vec[index]);
        }
        vec rhos(lLocIndices0.size());
        if(which_structfun == 0) {
/*            BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax); */
            BarnesStructure structure = BarnesStructure( dh, dz, dw);
            rhos = structure.corr_background(p1, p2);
        }
        else if(which_structfun == 1) {
/*            MixAStructure structure = MixAStructure( dh, dz, dw, hmax); */
            MixAStructure structure = MixAStructure( dh, dz, dw);
            rhos = structure.corr_background(p1, p2);
        } 
        for(int i = 0; i < lLocIndices0.size(); i++) {
            int index = lLocIndices0[i];
            if(gridpp::is_valid(pobs[index][0])) {
                if(rhos[i] > 0) {
                    lRhos0.push_back(std::pair<float,int>(rhos[i], i));
                }
            }
        }

        // Make sure we don't use too many observations
        arma::vec lRhos;
        if(max_points > 0 && lRhos0.size() > max_points) {
            // If sorting is enabled and we have too many locations, then only keep the best ones based on rho.
            // Otherwise, just use the last locations added
            lRhos = arma::vec(max_points);
            std::sort(lRhos0.begin(), lRhos0.end(), ::sort_pair_first<float,int>());
            for(int i = 0; i < max_points; i++) {
                // The best values start at the end of the array
                int index = lRhos0[lRhos0.size() - 1 - i].second;
                lLocIndices.push_back(lLocIndices0[index]);
                lRhos(i) = lRhos0[lRhos0.size() - 1 - i].first;
            }
        }
        else {
            lRhos = arma::vec(lRhos0.size());
            for(int i = 0; i < lRhos0.size(); i++) {
                int index = lRhos0[i].second;
                lLocIndices.push_back(lLocIndices0[index]);
                lRhos(i) = lRhos0[i].first;
            }
        }

        int lS = lLocIndices.size();
        if(lS == 0) {
            // If we have too few observations though, then use the background
            continue;
        }

        vectype lElevs(lS);
        vectype lLafs(lS);
        for(int i = 0; i < lLocIndices.size(); i++) {
            int index = lLocIndices[i];
            lElevs[i] = pelevs[index];
            lLafs[i] = plafs[index];
        }

        // Innovation: Observation - Background_r
        mattype lInnov(lS, nValidEns, arma::fill::zeros);
        // lR_dd: Observation error correlation matrix
        mattype lR_dd(lS, lS, arma::fill::zeros);
        // lCorr1D: background correlations between yth gridpoint and observations
        mattype lCorr1D(1, lS, arma::fill::zeros); 
        // lCorr2D: background correlations among observations
        mattype lCorr2D(lS, lS, arma::fill::zeros); 
        for(int i = 0; i < lS; i++) {
            // lR_dd: diagonal observation error correlation matrix
            lR_dd(i, i) = var_ratios_or;
            // lCorr1D between yth gridpoint and ith observation computed before
            lCorr1D(0, i) = lRhos(i);
            // compute lInnov
            int index = lLocIndices[i];
            for(int e = 0; e < nValidEns; e++) {
                int ei = validEns[e];
                lInnov(i, ei) = pobs[index][ei] - pbackground_r[index][ei];
            }
            // compute lCorr2D
            Point p1 = point_vec[index];
            std::vector<Point> p2(lS, Point(0, 0));
            for(int j = 0; j < lS; j++) {
                int index_j = lLocIndices[j];
                p2[j] = point_vec[index_j];
            }
            vec corr(lS);
            if(which_structfun == 0) {
/*                BarnesStructure structure = BarnesStructure( dh, dz, dw, hmax); */
                BarnesStructure structure = BarnesStructure( dh, dz, dw);
                corr = structure.corr_background(p1, p2);
            }
            else if(which_structfun == 1) {
/*                MixAStructure structure = MixAStructure( dh, dz, dw, hmax); */
                MixAStructure structure = MixAStructure( dh, dz, dw);
                corr = structure.corr_background(p1, p2);
            } 
            /* vec corr = structure.corr(p1, p2); */
            for(int j = 0; j < lS; j++) 
                lCorr2D(i, j) = corr[j];
        } // end loop over closer observations
        // lK(1, lS): Kalman gain
        mattype lK = lCorr1D * arma::inv(lCorr2D + lR_dd);
        // dx(1, nValidEns): analysis increment 
        mattype dx = std_ratios_lr * weight * (lK * lInnov);
        ///////////////////////////////
        // Anti-extrapolation filter //
        ///////////////////////////////
        if(!allow_extrapolation) {
            arma::rowvec maxIncEns = arma::max(lInnov);
            arma::rowvec minIncEns = arma::min(lInnov);

            for(int e = 0; e < nValidEns; e++) {
                float increment = dx[e];
                float maxInc = maxIncEns[e];
                float minInc = minIncEns[e];
                if(maxInc > 0 && increment > maxInc) {
                   increment = maxInc;
                }
                else if(maxInc < 0 && increment > 0) {
/*                   increment = maxInc; */
                   increment = 0; 
                }
                else if(minInc < 0 && increment < minInc) {
                   increment = minInc;
                }
                else if(minInc > 0 && increment < 0) {
/*                   increment = minInc; */
                   increment = 0;
                }
                dx[e] = increment;
            }
        }

        // Compute the analysis as the updated background_l
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            output[y][ei] = background_l[y][ei] + dx[e];
        }
    } // end loop over gridpoint 
    return output;
} // end of R_optimal_interpolation_ensi_staticcorr_lr
