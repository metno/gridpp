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

vec3 gridpp::optimal_interpolation_ensi_multi(const gridpp::Grid& bgrid,
        const vec2& bratios,
        const vec3& background, 
        const vec3& background_corr, 
        const gridpp::Points& points,
        const vec2& pobs,
        const vec& pratios, 
        const vec2& pbackground, 
        const vec2& pbackground_corr,
        const gridpp::StructureFunction& structure,
        int max_points,
        bool dynamic_correlations,
        bool allow_extrapolation) {
    double s_time = gridpp::clock();

    // Check input data
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");

    int nS = points.size();
    if(nS == 0)
        return background;

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
    int nE = background[0][0].size();
    if(background.size() != nY || background[0].size() != nX) {
        std::stringstream ss;
        ss << "Input left field (" << background.size() << "," << background[0].size() << "," << background[0][0].size() << ") is not the same size as the grid (" << nY << "," << nX << "," << nE << ")";
        throw std::invalid_argument(ss.str());
    }
    if(background_corr.size() != nY || background_corr[0].size() != nX || background_corr[0][0].size() != nE) {
        std::stringstream ss;
        ss << "Input LEFT field (" << background_corr.size() << "," << background_corr[0].size() << "," << background_corr[0][0].size() << ") is not the same size as the grid (" << nY << "," << nX << "," << nE << ")";
        throw std::invalid_argument(ss.str());
    }
    if(bratios.size() != nY || bratios[0].size() != nX) {
        std::stringstream ss;
        ss << "Input bratios field (" << bratios.size() << "," << bratios[0].size() << ") is not the same size as the grid (" << nY << "," << nX << ")";
        throw std::invalid_argument(ss.str());
    }
    if(pbackground.size() != nS || pbackground[0].size() != nE) {
        std::stringstream ss;
        ss << "Input right field at observation location (" << pbackground.size() << "," << pbackground[0].size() << ") and points (" << nS << "," << nE << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    if(pbackground_corr.size() != nS || pbackground_corr[0].size() != nE) {
        std::stringstream ss;
        ss << "Input RIGHT field at observation location (" << pbackground_corr.size() << "," << pbackground_corr[0].size() << ") and points (" << nS << "," << nE << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    if(pobs.size() != nS || pobs[0].size() != nE) {
        std::stringstream ss;
        ss << "Observations (" << pobs.size() << "," << pobs[0].size() << ") and points (" << nS << "," << nE << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    if(pratios.size() != nS) {
        std::stringstream ss;
        ss << "Ratios (" << pratios.size() << ") and points (" << nS << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
 
    gridpp::Points bpoints = bgrid.to_points();
    vec2 background1 = gridpp::init_vec2(nY * nX, nE);
    vec2 background_corr1 = gridpp::init_vec2(nY * nX, nE);
    vec bratios1(nY * nX);
    int count = 0;
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            bratios1[count] = bratios[y][x];
            for(int e = 0; e < nE; e++) {
                background1[count][e] = background[y][x][e];
                background_corr1[count][e] = background_corr[y][x][e];
            }
            count++;
        }
    }
    vec2 output1 = gridpp::init_vec2(nY * nX, nE); 
    if(dynamic_correlations) {
        output1 = optimal_interpolation_ensi_multi_ebe(bpoints, bratios1, background1, background_corr1, points, pobs, pratios, pbackground, pbackground_corr, structure, max_points, allow_extrapolation);
    }
    else {
        output1 = optimal_interpolation_ensi_multi_ebesc(bpoints, bratios1, background1, points, pobs, pratios, pbackground, structure, max_points, allow_extrapolation);
    }
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
}

// ensemble member by ensemble member (ebe)
vec2 gridpp::optimal_interpolation_ensi_multi_ebe(const gridpp::Points& bpoints,
        const vec& bratios,
        const vec2& background, 
        const vec2& background_corr, 
        const gridpp::Points& points,
        const vec2& pobs,
        const vec& pratios,
        const vec2& pbackground, 
        const vec2& pbackground_corr,
        const gridpp::StructureFunction& structure,
        int max_points,
        bool allow_extrapolation) {
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");
    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background and observations points must be of same coorindate type (lat/lon or x/y)");
    }
    if(background.size() != bpoints.size())
        throw std::invalid_argument("Input background field is not the same size as the grid");
    if(background_corr.size() != bpoints.size())
        throw std::invalid_argument("Input background_corr field is not the same size as the grid");
    if(bratios.size() != bpoints.size())
        throw std::invalid_argument("Bratios and grid size mismatch");
    if(pobs.size() != points.size())
        throw std::invalid_argument("Observations and points exception mismatch");
    if(pratios.size() != points.size())
        throw std::invalid_argument("Pratios and points size mismatch");
    if(pbackground.size() != points.size())
        throw std::invalid_argument("Background and points size mismatch");
    if(pbackground_corr.size() != points.size())
        throw std::invalid_argument("Background_corr and points size mismatch");

    float default_min_std = 0.0013;
    int nS = points.size();
    if(nS == 0)
        return background;

    int mY = -1;  // Write debug information for this station index
    bool diagnose = false;

    int nY = background.size();
    int nEns = background[0].size();

    // Prepare output matrix
    float missing_value = -99999.999;
/*    vec2 output = gridpp::init_vec2(nY, nEns, missing_value); */
    vec2 output = background;

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
            float value_l = background[y][e];
            float value_L = background_corr[y][e];
            if(!gridpp::is_valid(value_l) || !gridpp::is_valid(value_L))
                numInvalid++;
        }
        for(int i = 0; i < nS; i++) {
            float value_r = pbackground[i][e];
            float value_R = pbackground_corr[i][e];
            if(!gridpp::is_valid(value_r) || !gridpp::is_valid(value_R))
                numInvalid++;
        }
        if(numInvalid == 0) {
            validEns.push_back(e);
            nValidEns++;
        }
    }
    if(nValidEns == 0)
        return background;

    // gZ_R(nY, nValidEns): used to compute ensemble-based background correlations i) between yth gridpoint and observations ii) among observations
    vec2 gZ_R = gridpp::init_vec2(nY, nValidEns); // useful to compute dynamical correlations
    for(int i = 0; i < nS; i++) {
        vec pbackgroundValid_R(nValidEns);
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            pbackgroundValid_R[e] = pbackground_corr[i][ei];
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
        float std_ratios_lr = bratios[y];
        float lat = blats[y];
        float lon = blons[y];
        float elev = belevs[y];
        float laf = blafs[y];
        Point p1 = bpoints.get_point(y);
        float localizationRadius = structure.localization_distance(p1);
 
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
        vec rhos = structure.corr_background(p1, p2); 
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
            backgroundValid_L[e] = background_corr[y][ei];
        }
        float mean = gridpp::calc_statistic(backgroundValid_L, gridpp::Mean);
        float std = gridpp::calc_statistic(backgroundValid_L, gridpp::Std);
        if(gridpp::is_valid(mean) && gridpp::is_valid(std) && std > default_min_std) {
            for(int e = 0; e < nValidEns; e++) 
                lX_L(0,e) = 1 / sqrt(nValidEns-1) * (backgroundValid_L[e] - mean) / std;
        }
        // lZ_R: used to compute ensemble-based background correlations i) between yth gridpoint and observations ii) among observations
        mattype lZ_R(lS, nValidEns);
        // Innovation: Observation - background
        mattype lInnov(lS, nValidEns, arma::fill::zeros);
        // lR_dd: Observation error correlation matrix
        mattype lR_dd(lS, lS, arma::fill::zeros);
        // lLoc1D: localization for ensemble-based background correlations between yth gridpoint and observations
        mattype lLoc1D(1, lS, arma::fill::zeros); 
        // lLoc2D: localization for ensemble-based background correlations among observations
        mattype lLoc2D(lS, lS, arma::fill::zeros); 
        for(int i = 0; i < lS; i++) {
            int index = lLocIndices[i];
            // lR_dd: diagonal observation error correlation matrix
            lR_dd(i, i) = pratios[index];
            // lLoc1D between yth gridpoint and ith observation computed before
            lLoc1D(0, i) = lRhos(i);
            // compute lZ_R and lInnov
            for(int e = 0; e < nValidEns; e++) {
                lZ_R(i, e) = gZ_R[index][e];
                int ei = validEns[e];
                lInnov(i, ei) = pobs[index][ei] - pbackground[index][ei];
            }
            // compute lLoc2D
            Point p1 = point_vec[index];
            std::vector<Point> p2(lS, Point(0, 0));
            for(int j = 0; j < lS; j++) {
                int index_j = lLocIndices[j];
                p2[j] = point_vec[index_j];
            }
            vec corr = structure.corr(p1, p2);
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
        mattype dx = std_ratios_lr * (lK * lInnov);
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

        // Compute the analysis as the updated background
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            output[y][ei] = background[y][ei] + dx[e];
        }
        // debug
/*        for(int i = 0; i < lS; i++) {
            // compute lZ_R and lInnov
            int index = lLocIndices[i];
            std::cout << i << " backg_r obs " << pbackground[index][0] << " " << pobs[index][0] << std::endl;
        } // end loop over closer observations
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            std::cout << ei << " backg_l analysis " << background[y][ei] << " " << output[y][ei] << std::endl;
        } */

    } // end loop over gridpoint 
    return output;
} // end of optimal_interpolation_ensi_multi_ebe

// ensemble member by ensemble member with static correlations (ebesc)
vec2 gridpp::optimal_interpolation_ensi_multi_ebesc(const gridpp::Points& bpoints,
        const vec& bratios,
        const vec2& background, 
        const gridpp::Points& points,
        const vec2& pobs,  
        const vec& pratios,
        const vec2& pbackground, 
        const gridpp::StructureFunction& structure, 
        int max_points,
        bool allow_extrapolation) {
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");
    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background and observations points must be of same coorindate type (lat/lon or x/y)");
    }
    if(background.size() != bpoints.size())
        throw std::invalid_argument("Input background field is not the same size as the grid");
    if(bratios.size() != bpoints.size())
        throw std::invalid_argument("Bratios and grid size mismatch");
    if(pobs.size() != points.size())
        throw std::invalid_argument("Observations and points exception mismatch");
    if(pratios.size() != points.size())
        throw std::invalid_argument("Pratios and points size mismatch");
    if(pbackground.size() != points.size())
        throw std::invalid_argument("Background and points size mismatch");

    int nS = points.size();
    if(nS == 0)
        return background;

    int mY = -1;  // Write debug information for this station index
    bool diagnose = false;

    int nY = background.size();
    int nEns = background[0].size();

    // Prepare output matrix
    float missing_value = -99999.999;
    vec2 output = background;

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
            float value_l = background[y][e];
            if(!gridpp::is_valid(value_l))
                numInvalid++;
        }
        for(int i = 0; i < nS; i++) {
            float value_r = pbackground[i][e];
            if(!gridpp::is_valid(value_r))
                numInvalid++;
        }
        if(numInvalid == 0) {
            validEns.push_back(e);
            nValidEns++;
        }
    }
    if(nValidEns == 0)
        return background;

    // This causes segmentation fault when building the gridpp pypi package
    // 1) Tested removing num_condition_warning and num_real_part_warning from the parallel loop
    //    but it doesnt seem to help
    // #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float std_ratios_lr = bratios[y];
        float lat = blats[y];
        float lon = blons[y];
        float elev = belevs[y];
        float laf = blafs[y];
        Point p1 = bpoints.get_point(y);
        float localizationRadius = structure.localization_distance(p1);
 
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
        vec rhos = structure.corr_background(p1, p2);
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

        // Innovation: Observation - background
        mattype lInnov(lS, nValidEns, arma::fill::zeros);
        // lR_dd: Observation error correlation matrix
        mattype lR_dd(lS, lS, arma::fill::zeros);
        // lCorr1D: background correlations between yth gridpoint and observations
        mattype lCorr1D(1, lS, arma::fill::zeros); 
        // lCorr2D: background correlations among observations
        mattype lCorr2D(lS, lS, arma::fill::zeros); 
        for(int i = 0; i < lS; i++) {
            int index = lLocIndices[i];
            // lR_dd: diagonal observation error correlation matrix
            lR_dd(i, i) = pratios[index];
            // lCorr1D between yth gridpoint and ith observation computed before
            lCorr1D(0, i) = lRhos(i);
            // compute lInnov
            for(int e = 0; e < nValidEns; e++) {
                int ei = validEns[e];
                lInnov(i, ei) = pobs[index][ei] - pbackground[index][ei];
            }
            // compute lCorr2D
            Point p1 = point_vec[index];
            std::vector<Point> p2(lS, Point(0, 0));
            for(int j = 0; j < lS; j++) {
                int index_j = lLocIndices[j];
                p2[j] = point_vec[index_j];
            }
            vec corr = structure.corr(p1, p2);
            for(int j = 0; j < lS; j++) 
                lCorr2D(i, j) = corr[j];
        } // end loop over closer observations
        // lK(1, lS): Kalman gain
        mattype lK = lCorr1D * arma::inv(lCorr2D + lR_dd);
        // dx(1, nValidEns): analysis increment 
        mattype dx = std_ratios_lr * (lK * lInnov);
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

        // Compute the analysis as the updated background
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            output[y][ei] = background[y][ei] + dx[e];
        }
    } // end loop over gridpoint 
    return output;
} // end of optimal_interpolation_ensi_multi_ebesc

// use the ensemble mean (utem) 
vec2 gridpp::optimal_interpolation_ensi_multi_utem(const gridpp::Points& bpoints,
            const vec& bratios,
            const vec2& background,
            const vec2& background_corr,
            const gridpp::Points& points,
            const vec2& pobs,
            const vec& pratios,
            const vec2& pbackground,
            const vec2& pbackground_corr,
            const gridpp::StructureFunction& structure,
            int max_points,
            bool allow_extrapolation) {
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");
    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background and observations points must be of same coorindate type (lat/lon or x/y)");
    }
    if(background.size() != bpoints.size())
        throw std::invalid_argument("Input field is not the same size as the grid");
    if(background_corr.size() != bpoints.size())
        throw std::invalid_argument("Input background_corr field is not the same size as the grid");
    if(bratios.size() != bpoints.size())
        throw std::invalid_argument("Bratios and grid size mismatch");
    if(pobs.size() != points.size())
        throw std::invalid_argument("Observations and points exception mismatch");
    if(pratios.size() != points.size())
        throw std::invalid_argument("Pratios and points size mismatch");
    if(pbackground.size() != points.size())
        throw std::invalid_argument("Background and points size mismatch");
    if(pbackground_corr.size() != points.size())
        throw std::invalid_argument("Background_corr and points size mismatch");

    float default_min_std = 0.0013;
    int nS = points.size();
    if(nS == 0)
        return background;

    int mY = -1;  // Write debug information for this station index
    bool diagnose = false;

    int nY = background.size();
    int nEns = background[0].size();

    // Prepare output matrix
    vec2 output = background;

    int num_condition_warning = 0;
    int num_real_part_warning = 0;

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
            float value_l = background[y][e];
            float value_L = background_corr[y][e];
            if(!gridpp::is_valid(value_l) || !gridpp::is_valid(value_L))
                numInvalid++;
        }
        for(int i = 0; i < nS; i++) {
            float value_r = pbackground[i][e];
            float value_R = pbackground_corr[i][e];
            if(!gridpp::is_valid(value_r) || !gridpp::is_valid(value_R))
                numInvalid++;
        }
        if(numInvalid == 0) {
            validEns.push_back(e);
            nValidEns++;
        }
    }
    if(nValidEns == 0)
        return background;

    // Compute Y
    vec2 gY = gridpp::init_vec2(nY, nValidEns);
    vec2 gY_corr = gridpp::init_vec2(nY, nValidEns);
    vec gYhat(nS);
    float const_fact = 1 / sqrt(nValidEns-1);
    for(int i = 0; i < nS; i++) {
        vec pbackgroundValid(nValidEns);
        vec pbackgroundValid_corr(nValidEns);
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            pbackgroundValid[e] = pbackground[i][ei];
            pbackgroundValid_corr[e] = pbackground_corr[i][ei];
        }
        float mean = gridpp::calc_statistic(pbackgroundValid, gridpp::Mean);
        if(!gridpp::is_valid(mean)) {
            for(int e = 0; e < nValidEns; e++) 
                gY[i][e] = 0;
        }
        else {
            for(int e = 0; e < nValidEns; e++) 
                gY[i][e] = pbackgroundValid[e] - mean;
        }
        gYhat[i] = mean;
        float mean_corr = gridpp::calc_statistic(pbackgroundValid_corr, gridpp::Mean);
        float std_corr = gridpp::calc_statistic(pbackgroundValid_corr, gridpp::Std);
        if(!gridpp::is_valid(mean_corr) || !gridpp::is_valid(std_corr)) {
            for(int e = 0; e < nValidEns; e++) 
                gY_corr[i][e] = 0;
        }
        else {
            if(std <= default_min_std) {
                for(int e = 0; e < nValidEns; e++) 
                    gY_corr[i][e] = 0;
            }
            else {
                for(int e = 0; e < nValidEns; e++) 
                    gY_corr[i][e] = const_fact * (pbackgroundValid_corr[e] - mean_corr) / std_corr;
            }
        }
    }

    // This causes segmentation fault when building the gridpp pypi package
    // 1) Tested removing num_condition_warning and num_real_part_warning from the parallel loop
    //    but it doesnt seem to help
    // #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float std_ratios_lr = bratios[y];
        float lat = blats[y];
        float lon = blons[y];
        float elev = belevs[y];
        float laf = blafs[y];
        Point p1 = bpoints.get_point(y);
        float localizationRadius = structure.localization_distance(p1);

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
        vec rhos = structure.corr_background(p1, p2);
        for(int i = 0; i < lLocIndices0.size(); i++) {
            int index = lLocIndices0[i];
            if(gridpp::is_valid(pobs[index])) {
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

        vectype lObs(lS);
        vectype lElevs(lS);
        vectype lLafs(lS);
        for(int i = 0; i < lLocIndices.size(); i++) {
            int index = lLocIndices[i];
            lObs[i] = pobs[index];
            lElevs[i] = pelevs[index];
            lLafs[i] = plafs[index];
        }

        // Compute Y (model at obs-locations)
        mattype lY(lS, nValidEns);
        mattype lY_corr(lS, nValidEns);
        vectype lYhat(lS);

        for(int i = 0; i < lS; i++) {
            // Use the nearest neighbour for this location
            int index = lLocIndices[i];
            for(int e = 0; e < nValidEns; e++) {
                lY(i, e) = gY[index][e];
                lY_corr(i, e) = gY_corr[index][e];
            }
            lYhat[i] = gYhat[index];
        }

        // Compute Rinv
        mattype Rinv(lS, lS, arma::fill::zeros);
        for(int i = 0; i < lS; i++) {
            int index = lLocIndices[i];
            Rinv(i, i) = lRhos[i] / pratios[index];
        }

        // Compute C matrix
        // k x nS * nS x nS
        mattype C(nValidEns, lS);
        C = lY_corr.t() * Rinv;

        mattype Pinv(nValidEns, nValidEns);
//        float diag = 1 / delta * (nValidEns - 1);

        Pinv = C * lY_corr + 1 * arma::eye<mattype>(nValidEns, nValidEns);
        float cond = arma::rcond(Pinv);
        if(cond <= 0) {
            num_condition_warning++;
            continue;
        }

        // Compute sqrt of matrix. Armadillo 6.6 has this function, but on many systems this
        // is not available. Therefore, compute sqrt using the method found in 6.6
        // cxtype Wcx(nValidEns, nValidEns);
        // status = arma::sqrtmat(Wcx, (nValidEns - 1) * P);
        // mattype W = arma::real(Wcx);

        mattype P = arma::inv(Pinv);
        vectype eigval;
        mattype eigvec;
        bool status = arma::eig_sym(eigval, eigvec, (nValidEns - 1) * P);
        if(!status) {
            std::cout << "Cannot find eigenvector:" << std::endl;
            std::cout << "Lat: " << lat << std::endl;
            std::cout << "Lon: " << lon << std::endl;
            std::cout << "Elev: " << elev << std::endl;
            std::cout << "Laf: " << laf << std::endl;
            std::cout << "Pinv" << std::endl;
            print_matrix<mattype>(Pinv);
            std::cout << "P" << std::endl;
            print_matrix<mattype>(P);
            std::cout << "Y:" << std::endl;
            print_matrix<mattype>(lY_corr);
            std::cout << "lObs:" << std::endl;
            print_matrix<mattype>(lObs);
            std::cout << "Yhat" << std::endl;
            print_matrix<mattype>(lYhat);
        }
        eigval = sqrt(eigval);
        mattype Wcx = eigvec * arma::diagmat(eigval) * eigvec.t();
        mattype W = arma::real(Wcx);

        if(W.n_rows == 0) {
            num_real_part_warning++;
            continue;
        }

        // Compute PC
        mattype PC(nValidEns, lS);
        PC = P * C;

        // Compute w
        vectype w(nValidEns);
        if(diagnose)
            w = PC * (arma::ones<vectype>(lS));
        else
            w = PC * (lObs - lYhat);

        // Compute X (perturbations about model mean)
        vectype X(nValidEns);
        vectype X_corr(nValidEns);
        float total = 0;
        int count = 0;
        float total_corr = 0;
        int count_corr = 0;
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            float value = background[y][ei];
            if(gridpp::is_valid(value)) {
                X(e) = value;
                total += value;
                count++;
            }
            float value_corr = background_corr[y][ei];
            if(gridpp::is_valid(value_corr)) {
                X_corr(e) = value_corr;
                total_corr += value_corr;
                count_corr++;
            }
        }
        float ensMean = total / count;
        float ensStd = gridpp::calc_statistic(X, gridpp::Std);
        float ensMean_corr = total_corr / count_corr;
        float ensStd_corr = gridpp::calc_statistic(X_corr, gridpp::Std);
        for(int e = 0; e < nValidEns; e++) {
            X(e) -= ensMean;
            float value_corr = X_corr(e);
            if (ensStd > 0) {
                X_corr(e) = const_fact * (value_corr - ensMean_corr) / ensStd_corr;
            }
            else {
                X_corr(e) = 0;
            }
        }

        // Add w to W
        for(int e = 0; e < nValidEns; e++) {
            for(int e2 = 0; e2 < nValidEns; e2 ++) {
                W(e, e2) = ensStd * W(e, e2) + std_ratios_lr * w(e) ;
            }
        }

        // Write debugging information
        if(y == mY) {
            std::cout << "Lat: " << lat << std::endl;
            std::cout << "Lon: " << lon << " " << lat << " " << std::endl;
            std::cout << "Elev: " << elev << std::endl;
            std::cout << "Laf: " << laf << std::endl;
            std::cout << "Num obs: " << lS << std::endl;
            std::cout << "Num ens: " << nValidEns << std::endl;
            std::cout << "rhos" << std::endl;
            print_matrix<mattype>(lRhos);
            std::cout << "P" << std::endl;
            print_matrix<mattype>(P);
            std::cout << "C" << std::endl;
            print_matrix<mattype>(C);
            std::cout << "C * lY_corr" << std::endl;
            print_matrix<mattype>(C * lY_corr);
            std::cout << "PC" << std::endl;
            print_matrix<mattype>(PC);
            std::cout << "W" << std::endl;
            print_matrix<mattype>(W);
            std::cout << "w" << std::endl;
            print_matrix<mattype>(w);
            std::cout << "Y:" << std::endl;
            print_matrix<mattype>(lY_corr);
            std::cout << "Yhat" << std::endl;
            print_matrix<mattype>(lYhat);
            std::cout << "lObs" << std::endl;
            print_matrix<mattype>(lObs);
            std::cout << "lObs - Yhat" << std::endl;
            print_matrix<mattype>(lObs - lYhat);
            std::cout << "X" << std::endl;
            print_matrix<mattype>(X);
            std::cout << "elevs" << std::endl;
            print_matrix<mattype>(lElevs);
            std::cout << "lafs" << std::endl;
            print_matrix<mattype>(lLafs);
            std::cout << "Analysis increment:" << std::endl;
            print_matrix<mattype>(X.t() * W);
            std::cout << "My: " << arma::mean(arma::dot(lObs - lYhat, lRhos) / lS) << std::endl;
        }

        // Compute analysis
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            float total = 0;
            for(int k = 0; k < nValidEns; k++) {
                total += X_corr(k) * W(k, e);
            }

            float currIncrement = total;

            float raw = ensMean;

            ///////////////////////////////
            // Anti-extrapolation filter //
            ///////////////////////////////
            if(!allow_extrapolation) {
                // Don't allow a final increment that is larger than any increment
                // at station points
                float maxInc = arma::max(lObs - (lY[e] + lYhat));
                float minInc = arma::min(lObs - (lY[e] + lYhat));
                if(y == mY) {
                    std::cout << "Increments: " << maxInc << " " << minInc << " " << currIncrement << std::endl;
                }

                // The increment for this member. currIncrement is the increment relative to
                // ensemble mean
                float memberIncrement = currIncrement - X(e);
                // Adjust increment if it gives a member increment that is outside the range
                // of the observation increments
                if(y == mY) {
                    std::cout << "Analysis increment: " << memberIncrement << " " << ensMean << " " << currIncrement << " " << X(e) << std::endl;
                }
                if(maxInc > 0 && memberIncrement > maxInc) {
                    currIncrement = maxInc + X(e);
                }
                else if(maxInc < 0 && memberIncrement > 0) {
                    currIncrement = 0 + X(e);
                }
                else if(minInc < 0 && memberIncrement < minInc) {
                    currIncrement = minInc + X(e);
                }
                else if(minInc > 0 && memberIncrement < 0) {
                    currIncrement = 0 + X(e);
                }
                if(y == mY) {
                    std::cout << "Final increment: " << currIncrement << " " << currIncrement - X(e) << std::endl;
                }
            }
            output[y][ei] = ensMean + currIncrement;
        }
    }

    if(num_condition_warning > 0) {
        std::stringstream ss;
        ss << "Condition number error in " << num_condition_warning << " points. Using raw values in those points.";
        gridpp::warning(ss.str());
    }
    if(num_real_part_warning > 0) {
        std::stringstream ss;
        ss << "Could not find the real part of W in " << num_real_part_warning << " points. Using raw values in those points.";
        gridpp::warning(ss.str());
    }
    return output;
}  // end of optimal_interpolation_ensi_multi_utem
