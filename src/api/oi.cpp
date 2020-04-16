#include "gridpp.h"
#include <math.h>
#include <algorithm>
#include <armadillo>
#include <assert.h>
#include <exception>

namespace {
    typedef arma::mat mattype;
    typedef arma::vec vectype;
    typedef arma::cx_mat cxtype;

    float calcRho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wlength);
    void check_vec(vec2 input, int Y, int X);
    void check_vec(vec input, int S);

    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };

    vec compute_background(const vec2& input, const gridpp::Grid& grid, const gridpp::Points& points, float elev_gradient);
}

vec2 gridpp::optimal_interpolation(const gridpp::Grid& bgrid,
        const vec2& input,
        const gridpp::Points& points,
        const vec& pobs,  // gObs
        const vec& pci,   // gCi
        const gridpp::StructureFunction& structure,
        int max_points,
        float elev_gradient,
        float epsilon) {
    double s_time = gridpp::util::clock();

    // Check input data
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");
    if(input.size() != bgrid.size()[0] || input[0].size() != bgrid.size()[1])
        throw std::runtime_error("Input field is not the same size as the grid");
    if(pobs.size() != points.size())
        throw std::runtime_error("Observations and points exception mismatch");
    if(pci.size() != points.size())
        throw std::runtime_error("Ci and points size mismatch");

    int nY = input.size();
    int nX = input[0].size();

    // Prepare output matrix
    vec2 output = gridpp::util::init_vec2(nY, nX);

    vec2 blats = bgrid.get_lats();
    vec2 blons = bgrid.get_lons();
    vec2 belevs = bgrid.get_elevs();
    vec2 blafs = bgrid.get_lafs();

    ////////////////////////////////////
    // Remove stations outside domain //
    ////////////////////////////////////
    ivec indices = points.get_in_domain_indices(bgrid);
    gridpp::Points points0;
    points0 = points.get_in_domain(bgrid);

    vec plats = points0.get_lats();
    vec plons = points0.get_lons();
    vec pelevs = points0.get_elevs();
    vec plafs = points0.get_lafs();
    int nS = plats.size();
    assert(indices.size() == nS);
    vec pobs0(nS);
    vec pci0(nS);
    for(int s = 0; s < nS; s++) {
        pobs0[s] = pobs[indices[s]];
        pci0[s] = pci[indices[s]];
    }

    std::stringstream ss;
    ss << "Number of observations: " << nS << std::endl;
    ss << "Number of gridpoints: " << nY << " " << nX;
    gridpp::util::debug(ss.str());

    float localizationRadius = structure.localization_distance();

    // Store the nearst gridpoint indicies for each observation
    std::vector<float> pYi(nS, gridpp::MV);
    std::vector<float> pXi(nS, gridpp::MV);
    for(int i = 0; i < nS; i++) {
        const std::vector<int> indices = bgrid.get_nearest_neighbour(plats[i], plons[i]);
        pYi[i] = indices[0];
        pXi[i] = indices[1];
    }

    // Transform the background
    // #pragma omp parallel for
    // for(int x = 0; x < nX; x++) {
    //     for(int y = 0; y < nY; y++) {
    //         float value = input[y][x];
    //         if(gridpp::util::is_valid(value))
    //             input[y][x] = transform.forward(value);
    //     }
    // }

    // Compute the background value at observation points (Y)
    vec gY = ::compute_background(input, bgrid, points, elev_gradient);

    #pragma omp parallel for
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            float lat = blats[y][x];
            float lon = blons[y][x];
            float elev = belevs[y][x];
            float laf = blafs[y][x];

            // Find observations within localization radius
            // TODO: Check that the chosen ones have elevation
            ivec lLocIndices0 = points0.get_neighbours(lat, lon, localizationRadius);
            if(lLocIndices0.size() == 0) {
                // If we have too few observations though, then use the background
                output[y][x] = input[y][x];
                continue;
            }
            std::vector<int> lLocIndices;
            lLocIndices.reserve(lLocIndices0.size());
            std::vector<std::pair<float,int> > lRhos0;

            // Calculate gridpoint to observation rhos
            lRhos0.reserve(lLocIndices0.size());
            for(int i = 0; i < lLocIndices0.size(); i++) {
                int index = lLocIndices0[i];
                Point p1(plats[index], plons[index], pelevs[index], plafs[index]);
                Point p2(lat, lon, elev, laf);
                float rho = structure.corr(p1, p2);
                if(rho > 0) {
                    lRhos0.push_back(std::pair<float,int>(rho, i));
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
            if(lS == 0) {
                // If we have too few observations though, then use the background
                output[y][x] = input[y][x];
                continue;
            }

            vectype lObs(lS);
            // Compute Y (model at obs-locations)
            vectype lY(lS);
            // Current grid-point to station error covariance matrix
            mattype lG(1, lS, arma::fill::zeros);
            // Station to station error covariance matrix
            mattype lP(lS, lS, arma::fill::zeros);
            // Station variance
            mattype lR(lS, lS, arma::fill::zeros);
            for(int i = 0; i < lS; i++) {
                int index = lLocIndices[i];
                lObs(i) = pobs0[index];
                lY(i) = gY[index];
                lR(i, i) = pci0[index];
                lG(0, i) = lRhos(i);
                Point p1(plats[index], plons[index], pelevs[index], plafs[index]);
                for(int j = 0; j < lS; j++) {
                    int index_j = lLocIndices[j];
                    Point p2(plats[index_j], plons[index_j], pelevs[index_j], plafs[index_j]);
                    lP(i, j) = structure.corr(p1, p2);
                }
            }
            mattype lGSR = lG * arma::inv(lP + epsilon * epsilon * lR);
            vectype dx = lGSR * (lObs - lY);
            output[y][x] = input[y][x] + dx[0];
        }
    }

    // Back-transform
    // #pragma omp parallel for
    // for(int x = 0; x < nX; x++) {
    //     for(int y = 0; y < nY; y++) {
    //         float value = (*output)(y, x, e);
    //         if(gridpp::util::is_valid(value))
    //             (*output)(y, x, e) = transform.backward(value);
    //     }
    // }
    std::cout << "OI total time: " << gridpp::util::clock() - s_time << std::endl;
    return output;
}

namespace {
    float calcRho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wlength) {
       float h = (iHDist/hlength);
       float rho = exp(-0.5 * h * h);
       if(gridpp::util::is_valid(vlength) && vlength > 0) {
          if(!gridpp::util::is_valid(iVDist)) {
             rho = 0;
          }
          else {
             float v = (iVDist/vlength);
             float factor = exp(-0.5 * v * v);
             rho *= factor;
          }
       }
       if(gridpp::util::is_valid(wlength) && wlength > 0) {
          float factor = exp(-0.5 * iLDist * iLDist / (wlength * wlength));
          rho *= factor;
       }
       return rho;
    }
    // Set up convenient functions for debugging in gdb
    // template<class Matrix>
    // void print_matrix(Matrix matrix) {
    //     matrix.print(std::cout);
    // }

    // template void print_matrix<CalibratorOi::mattype>(CalibratorOi::mattype matrix);
    // template void print_matrix<CalibratorOi::cxtype>(CalibratorOi::cxtype matrix);
    void check_vec(vec2 input, int Y, int X) {
        std::cout << input.size() << " " << input[0].size() << " " << Y << " " << X << std::endl;
        assert(input.size() == Y);
        for(int i = 0; i < input.size(); i++) {
            assert(input[i].size() == X);
            for(int j = 0; j < input[i].size(); j++) {
                assert(gridpp::util::is_valid(input[i][j]));
            }
        }
    }
    void check_vec(vec input, int S) {
        assert(input.size() == S);
        for(int i = 0; i < input.size(); i++) {
            assert(gridpp::util::is_valid(input[i]));
        }
    }
    vec compute_background(const vec2& input, const gridpp::Grid& grid, const gridpp::Points& points, float elev_gradient) {
        vec output(points.size());
        vec lats = points.get_lats();
        vec lons = points.get_lons();
        vec elevs = points.get_elevs();
        vec2 gelevs = grid.get_elevs();
        for(int i = 0; i < points.size(); i++) {
            ivec indices = grid.get_nearest_neighbour(lats[i], lons[i]);
            int y = indices[0];
            int x = indices[1];
            output[i] = input[y][x];
            if(gridpp::util::is_valid(elev_gradient) && elev_gradient != 0) {
                float nnElev = gelevs[y][x];
                assert(gridpp::util::is_valid(nnElev));
                assert(gridpp::util::is_valid(elevs[i]));
                float elevDiff = elevs[i] - nnElev;
                float elevCorr = elev_gradient * elevDiff;
                output[i] += elevCorr;
            }
        }
        return output;
    }
}
/*
float transform(float iValue) const {
   if(mTransformType == TransformTypeNone)
      return iValue;

   if(iValue <= 0)
      iValue = 0;
   if(mLambda == 0)
      return log(iValue);
   else
      return (pow(iValue, mLambda) - 1) / mLambda;
}

float CalibratorOi::invTransform(float iValue) const {
   if(mTransformType == TransformTypeNone)
      return iValue;

   float rValue = 0;
   if(mLambda == 0)
      rValue = exp(iValue);
   else {
      if(iValue < -1.0 / mLambda) {
         iValue = -1.0 / mLambda;
      }
      rValue = pow(1 + mLambda * iValue, 1 / mLambda);
   }
   if(rValue <= 0)
      rValue = 0;
   return rValue;
}
*/

