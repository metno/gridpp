#include "gridpp.h"
#include <math.h>
#include <algorithm>
#include <armadillo>
#include <assert.h>
#include <exception>
#include <boost/math/distributions/normal.hpp>

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

vec2 gridpp::staticcorr_points(const gridpp::Points& points,
        const gridpp::Points& knots,
        const gridpp::StructureFunction& structure,
        int max_points) {

    double s_time = gridpp::clock();

    // Check input data
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");

    if(points.get_coordinate_type() != knots.get_coordinate_type()) {
        throw std::invalid_argument("Both background grid and observations points must be of same coordinate type (lat/lon or x/y)");
    }

    int nY = points.size();
    int nS = knots.size();

    vec2 output = gridpp::init_vec2(nY, nS, 0);

    vec blats = points.get_lats();
    vec blons = points.get_lons();
    vec belevs = points.get_elevs();
    vec blafs = points.get_lafs();

    vec plats = knots.get_lats();
    vec plons = knots.get_lons();
    vec pelevs = knots.get_elevs();
    vec plafs = knots.get_lafs();

    // Create all objects of type Point (to save time on creation later)
    std::vector<Point> point_vec;
    point_vec.reserve(nS);
    for(int i = 0; i < nS; i++) {
        point_vec.push_back(knots.get_point(i));
    }

    #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float lat = blats[y];
        float lon = blons[y];
        Point p1 = points.get_point(y);
        float localizationRadius = structure.localization_distance(p1);

        // Find observations within localization radius
        // TODO: Check that the chosen ones have elevation
        ivec lLocIndices0 = knots.get_neighbours(lat, lon, localizationRadius);
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
            if(rhos[i] > 0) {
                lRhos0.push_back(std::pair<float,int>(rhos[i], i));
            }
        }

        // Make sure we don't use too many observations
        arma::vec lRhos;
        if(max_points > 0 && lRhos0.size() > max_points) {
            // If we have too many locations, then only keep the best ones based on rho.
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
        if(lS > 0) {
            for(int i = 0; i < lS; i++) {
                int index = lLocIndices[i];
                output[y][index] = lRhos(i);
            }
        }
    } // End loop over points
    double e_time = gridpp::clock();
    // std::cout << count_stat << " " << e_time - s_time << " s" << std::endl;
    return output;
} // End function corr_points
