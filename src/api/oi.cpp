#include "gridpp.h"
#include <math.h>
#include <algorithm>
#include <armadillo>
#include <assert.h>

namespace {
    typedef arma::mat mattype;
    typedef arma::vec vectype;
    typedef arma::cx_mat cxtype;

    float calcRho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wmin);
    void check_vec(vec2 input, int Y, int X);
    void check_vec(vec input, int S);

    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };
}

int gridpp::optimal_interpolation(const vec2& input,
        const gridpp::Grid& bgrid,
        const vec& pobs,  // gObs
        const vec& pci,   // gCi
        const gridpp::Points& points,
        float minRho,
        float hlength,
        float vlength,
        float wmin,
        int maxPoints,
        float elevGradient,
        float epsilon,
        vec2& output) {
    double s_time = gridpp::util::clock();

    int nY = input.size();
    int nX = input[0].size();
    int nS = points.size();

    // Check input data
    check_vec(pobs, nS);
    check_vec(pci, nS);
    vec2 blats = bgrid.get_lats();
    vec2 blons = bgrid.get_lons();
    vec2 belevs = bgrid.get_elevs();
    vec2 blafs = bgrid.get_lafs();
    vec plats = points.get_lats();
    vec plons = points.get_lons();
    vec pelevs = points.get_elevs();
    vec plafs = points.get_lafs();

    // Prepare output matrix
    output.resize(nY);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nX);
    }

    // Estimate the grid spacing
    float gridSize = gridpp::KDTree::calc_distance(blats[0][0], blons[0][0], blats[1][0], blons[1][0]);
    std::stringstream ss;
    ss << "Estimated grid size: " << gridSize << " m" << std::endl;
    ss << "Number of observations: " << nS << std::endl;
    ss << "Number of gridpoints: " << nY << " " << nX;
    gridpp::util::debug(ss.str());

    // Loop over each observation, find the nearest gridpoint and place the obs into all gridpoints
    // in the vicinity of the nearest neighbour. This is only meant to be an approximation, but saves
    // considerable time instead of doing a loop over each grid point and each observation.

    // Store the indicies (into the gPoints array) that a gridpoint has available
    std::vector<std::vector<std::vector<int> > > gLocIndices; // Y, X, obs indices
    std::vector<float> pYi(nS, gridpp::MV);
    std::vector<float> pXi(nS, gridpp::MV);
    gLocIndices.resize(nY);
    for(int y = 0; y < nY; y++) {
        gLocIndices[y].resize(nX);
    }

    // Calculate the factor that the horizontal decorrelation scale should be multiplied by
    // to get the localization radius. For minRho=0.0013, the factor is 3.64
    float radiusFactor = sqrt(-2*log(minRho));

    // Spread each observation out to this many gridpoints from the nearest neighbour
    int gridpointRadius = radiusFactor * hlength / gridSize;

    // Vectorize the matrix of grid points
    vec blats0(nX * nY);
    vec blons0(nX * nY);
    int count = 0;
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            blats0[count] = blats[y][x];
            blons0[count] = blons[y][x];
            count++;
        }
    }

    // For each gridpoint, find which observations are relevant. Parse the observations and only keep
    // those that pass certain checks
    for(int i = 0; i < nS; i++) {
        const std::vector<int> indices = bgrid.get_nearest_neighbour(plats[i], plons[i]);
        int X = indices[1];
        int Y = indices[0];
        pYi[i] = Y;
        pXi[i] = X;

        // Check if the elevation of the station roughly matches the reference grid elevation
        bool hasValidElev = gridpp::util::is_valid(pelevs[i]);
        if(hasValidElev) {
            for(int y = std::max(0, Y - gridpointRadius); y < std::min(nY, Y + gridpointRadius); y++) {
                for(int x = std::max(0, X - gridpointRadius); x < std::min(nX, X + gridpointRadius); x++) {
                    gLocIndices[y][x].push_back(i);
                }
            }
        }
    }
    std::cout << "Done assigning observations to gridpoints " << gridpp::util::clock() - s_time << std::endl;

    // Transform the background
    // #pragma omp parallel for
    // for(int x = 0; x < nX; x++) {
    //     for(int y = 0; y < nY; y++) {
    //         float value = input[y][x];
    //         if(gridpp::util::is_valid(value))
    //             input[y][x] = gridpp::transform(value);
    //     }
    // }

    // Compute Y
    vec gY(nS);
    for(int i = 0; i < nS; i++) {
        float elevCorr = 0;
        if(gridpp::util::is_valid(elevGradient) && elevGradient != 0) {
            float nnElev = belevs[pYi[i]][pXi[i]];
            assert(gridpp::util::is_valid(nnElev));
            assert(gridpp::util::is_valid(pelevs[i]));
            float elevDiff = pelevs[i] - nnElev;
            elevCorr = elevGradient * elevDiff;
        }
        float value = input[pYi[i]][pXi[i]];
        // if(gridpp::util::is_valid(value)) {
        value += elevCorr;
        gY[i] = value;
    }

    #pragma omp parallel for
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            float lat = blats[y][x];
            float lon = blons[y][x];
            float elev = belevs[y][x];
            float laf = blafs[y][x];

            // Create list of locations for this gridpoint
            std::vector<int> lLocIndices0 = gLocIndices[y][x];
            if(lLocIndices0.size() == 0) {
                // If we have too few observations though, then use the background
                output[y][x] = input[y][x];
                continue;
            }
            std::vector<int> lLocIndices;
            lLocIndices.reserve(lLocIndices0.size());
            std::vector<std::pair<float,int> > lRhos0;
            lRhos0.reserve(lLocIndices0.size());
            for(int i = 0; i < lLocIndices0.size(); i++) {
                int index = lLocIndices0[i];
                float hdist = gridpp::KDTree::calc_distance(plats[index], plons[index], lat, lon);
                float vdist = gridpp::MV;
                if(gridpp::util::is_valid(pelevs[index] && gridpp::util::is_valid(elev)))
                    vdist = pelevs[index] - elev;
                float lafdist = 0;
                if(gridpp::util::is_valid(plafs[index]) && gridpp::util::is_valid(laf))
                    lafdist = plafs[index] - laf;
                float rho = ::calcRho(hdist, vdist, lafdist, hlength, vlength, wmin);
                int X = pXi[index];
                int Y = pYi[index];
                // Only include observations that are within the domain
                if(X > 0 && X < blats[0].size()-1 && Y > 0 && Y < blats.size()-1) {
                    if(rho > minRho) {
                        lRhos0.push_back(std::pair<float,int>(rho, i));
                    }
                }
            }

            arma::vec lRhos;
            if(lRhos0.size() > maxPoints) {
                // If sorting is enabled and we have too many locations, then only keep the best ones based on rho.
                // Otherwise, just use the last locations added
                lRhos = arma::vec(maxPoints);
                std::sort(lRhos0.begin(), lRhos0.end(), ::sort_pair_first<float,int>());
                for(int i = 0; i < maxPoints; i++) {
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
            vectype lElevs(lS);
            vectype lLafs(lS);
            for(int i = 0; i < lLocIndices.size(); i++) {
                int index = lLocIndices[i];
                lObs[i] = pobs[index];
                lElevs[i] = pelevs[index];
                lLafs[i] = plafs[index];
            }

            // Compute Y (model at obs-locations)
            vectype lY(lS);

            for(int i = 0; i < lS; i++) {
                // Use the nearest neighbour for this location
                int index = lLocIndices[i];
                lY(i) = gY[index];
            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // Single-member mode:                                                                //
            // Revert to static structure function when there is not enough ensemble information  //
            ////////////////////////////////////////////////////////////////////////////////////////
            // Current grid-point to station error covariance matrix
            mattype lG(1, lS, arma::fill::zeros);
            // Station to station error covariance matrix
            mattype lP(lS, lS, arma::fill::zeros);
            // Station variance
            mattype lR(lS, lS, arma::fill::zeros);
            for(int i = 0; i < lS; i++) {
                int index = lLocIndices[i];
                lR(i, i) = pci[index];
                float hdist = gridpp::KDTree::calc_distance(plats[index], plons[index], lat, lon);
                float vdist = gridpp::MV;
                if(gridpp::util::is_valid(pelevs[index] && gridpp::util::is_valid(elev)))
                    vdist = pelevs[index] - elev;
                float lafdist = 0;
                if(gridpp::util::is_valid(plafs[index]) && gridpp::util::is_valid(laf))
                    lafdist = plafs[index] - laf;
                float rho = ::calcRho(hdist, vdist, lafdist, hlength, vlength, wmin);
                lG(0, i) = rho;
                for(int j = 0; j < lS; j++) {
                    int index_j = lLocIndices[j];
                    float hdist = gridpp::KDTree::calc_distance(plats[index], plons[index], plats[index_j], plons[index_j]);
                    float vdist = gridpp::MV;
                    if(gridpp::util::is_valid(pelevs[index] && gridpp::util::is_valid(pelevs[index_j])))
                        vdist = pelevs[index] - pelevs[index_j];
                    float lafdist = 0;
                    if(gridpp::util::is_valid(plafs[index]) && gridpp::util::is_valid(laf))
                        lafdist = plafs[index] - plafs[index_j];

                    lP(i, j) = ::calcRho(hdist, vdist, lafdist, hlength, vlength, wmin);
                }
            }
            mattype lGSR;
            lGSR = lG * arma::inv(lP + epsilon * epsilon * lR);

            vectype currFcst = lY;
            vectype dx = lGSR * (lObs - currFcst);
            output[y][x] = input[y][x] + dx[0];
        }
    }

    // Back-transform
    /*
    #pragma omp parallel for
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            float value = (*output)(y, x, e);
            if(gridpp::util::is_valid(value))
                (*output)(y, x, e) = invTransform(value);
        }
    }
    */
    std::cout << "OI total time: " << gridpp::util::clock() - s_time << std::endl;
    return 0;

    // return optimal_interpolation_single_member(input, bgrid.get_lats(), bgrid.get_lons(), bgrid.get_elevs(), bgrid.get_lafs(), pobs, pci, points.get_lats(), points.get_lons(), points.get_elevs(), points.get_lafs(), minRho, hlength, vlength, wmin, maxPoints, elevGradient, epsilon, output);
}
/*
int gridpp::optimal_interpolation_single_member(const vec2& input,
        const vec2& blats,
        const vec2& blons,
        const vec2& belevs,
        const vec2& blafs,
        const vec& pobs,  // gObs
        const vec& pci,   // gCi
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
        vec2& output) {

    gridpp::Grid bgrid(blats, blons, belevs, blafs);
    gridpp::Points points(plats, plons, pelevs, plafs);
    return optimal_interpolation_single_member(input, bgrid, pobs, pci, points, minRho, hlength, vlength, wmin, maxPoints, elevGradient, epsilon, output);
}
*/
namespace {
    float calcRho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wmin) {
       float h = (iHDist/hlength);
       float rho = exp(-0.5 * h * h);
       if(gridpp::util::is_valid(vlength)) {
          if(!gridpp::util::is_valid(iVDist)) {
             rho = 0;
          }
          else {
             float v = (iVDist/vlength);
             rho *= exp(-0.5 * v * v);
          }
       }
       if(gridpp::util::is_valid(wmin)) {
          rho *= 1 - (1 - wmin) * std::abs(iLDist);
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

int gridpp::optimal_interpolation_ens(const vec3& input,
            const gridpp::Grid& bgrid,
            const vec& pobs,
            const vec& pci,
            const gridpp::Points& points,
            vec2& output) {

}
