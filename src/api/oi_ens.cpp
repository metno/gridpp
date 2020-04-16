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
}

template<class Matrix>
void print_matrix(Matrix matrix) {
       matrix.print(std::cout);
}

template void print_matrix< ::mattype>(::mattype matrix);
template void print_matrix< ::cxtype>(::cxtype matrix);

vec3 gridpp::optimal_interpolation_ens(const gridpp::Grid& bgrid,
        const vec3& input,
        const gridpp::Points& points,
        const vec& pobs,  // gObs
        const vec& pci,   // pci
        float min_rho,
        float hlength,
        float vlength,
        float wlength,
        int max_points,
        float elev_gradient,
        float epsilon) {
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");
    if(min_rho <= 0)
        throw std::invalid_argument("min_rho must be > 0");
    std::cout << bgrid.size()[0] << " " << bgrid.size()[1] << " " << input.size() << " " << input[0].size() << std::endl;
    if(input.size() != bgrid.size()[0] || input[0].size() != bgrid.size()[1])
        throw std::runtime_error("Input field is not the same size as the grid");
    if(pobs.size() != points.size())
        throw std::runtime_error("Observations and points exception mismatch");
    if(pci.size() != points.size())
        throw std::runtime_error("Ci and points size mismatch");

    double s_time = gridpp::util::clock();
    int mX = 0;
    int mY = 0;
    int mMinValidEns = 5;
    int numParameters = 2;
    float sigma = 0.5;
    float sigmac = 0.5;
    float delta = 1;
    bool mExtrapolate = false;
    bool mDiagnose = false;

    int nY = input.size();
    int nX = input[0].size();
    int nEns = input[0][0].size();
    int nS = points.size();

    // Check input data
    vec2 blats = bgrid.get_lats();
    vec2 blons = bgrid.get_lons();
    vec2 belevs = bgrid.get_elevs();
    vec2 blafs = bgrid.get_lafs();
    vec plats = points.get_lats();
    vec plons = points.get_lons();
    vec pelevs = points.get_elevs();
    vec plafs = points.get_lafs();

    // Prepare output matrix
    vec3 output;
    output.resize(nY);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nX);
        for(int x = 0; x < nX; x++) {
            output[y][x].resize(nEns);
        }
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
    // to get the localization radius. For min_rho=0.0013, the factor is 3.64
    float radiusFactor = sqrt(-2*log(min_rho));

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
    vec2 gY(nS);
    vec gYhat(nS);
    for(int i = 0; i < nS; i++) {
        gY[i].resize(nEns, 0);
        float elevCorr = 0;
        if(gridpp::util::is_valid(elev_gradient) && elev_gradient != 0) {
            float nnElev = belevs[pYi[i]][pXi[i]];
            assert(gridpp::util::is_valid(nnElev));
            assert(gridpp::util::is_valid(pelevs[i]));
            float elevDiff = pelevs[i] - nnElev;
            elevCorr = elev_gradient * elevDiff;
        }

        float total = 0;
        int count = 0;
        for(int e = 0; e < nEns; e++) {
            float value = input[pYi[i]][pXi[i]][e];
            if(gridpp::util::is_valid(value)) {
                value += elevCorr;
                gY[i][e] = value;
                total += value;
                count++;
            }
        }
        float mean = gridpp::MV;
        if(count > 0)
            mean = total / count;
        for(int e = 0; e < nEns; e++) {
            float value = gY[i][e];
            if(gridpp::util::is_valid(value) && gridpp::util::is_valid(mean)) {
                gY[i][e] -= mean;
            }
        }
        gYhat[i] = mean;
    }

    // Calculate number of valid members
    int nValidEns = 0;
    std::vector<int> validEns;
    for(int e = 0; e < nEns; e++) {
        int numInvalid = 0;
        for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
                float value = input[y][x][e];
                if(!gridpp::util::is_valid(value))
                    numInvalid++;
            }
        }
        if(numInvalid == 0) {
            validEns.push_back(e);
            nValidEns++;
        }
    }
    std::cout << "Number of valid ensemble members: " << nValidEns << std::endl;
    // TODO: Deal with single member mode
    bool singleMemberMode = nValidEns < mMinValidEns;

    // #pragma omp parallel for
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
                float rho = ::calcRho(hdist, vdist, lafdist, hlength, vlength, wlength);
                int X = pXi[index];
                int Y = pYi[index];
                // Only include observations that are within the domain
                if(X > 0 && X < blats[0].size()-1 && Y > 0 && Y < blats.size()-1) {
                    if(rho > min_rho) {
                        lRhos0.push_back(std::pair<float,int>(rho, i));
                    }
                }
            }

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
            mattype lY(lS, nValidEns);
            vectype lYhat(lS);

            for(int i = 0; i < lS; i++) {
                // Use the nearest neighbour for this location
                int index = lLocIndices[i];
                for(int e = 0; e < nValidEns; e++) {
                    int ei = validEns[e];
                    lY(i, e) = gY[index][ei];
                }
                lYhat[i] = gYhat[index];
            }

            // Compute Rinv
            mattype Rinv(lS, lS, arma::fill::zeros);
            if(numParameters == 2) {
                for(int i = 0; i < lS; i++) {
                    int index = lLocIndices[i];
                    Rinv(i, i) = lRhos[i] / (sigma * sigma * pci[index]);
                }
            }
            else if(numParameters == 3) {
                /*
                // Inverting the matrix is more complicated, since the radar observations
                // have covariances. Therefore invert the covariance matrix for the radar part and
                // insert the values into the big inverse matrix.
                // std::cout << "Computing R matrix" << std::endl;
                // R = get_precipitation_r(gRadarL, pci, lLocIndices, lRhos);
                // Compute little R
                std::vector<int> gRadarIndices;
                gRadarIndices.reserve(lS);
                std::vector<int> lRadarIndices;
                lRadarIndices.reserve(lS);
                for(int i = 0; i < lS; i++) {
                    int index = lLocIndices[i];
                    if(gRadarL[index] > 0) {
                        gRadarIndices.push_back(index);
                        lRadarIndices.push_back(i);
                    }
                }
                int lNumRadar = gRadarIndices.size();

                // Compute R tilde r
                mattype radarR(lNumRadar, lNumRadar, arma::fill::zeros);
                for(int i = 0; i < lNumRadar; i++) {
                    for(int j = 0; j < lNumRadar; j++) {
                        int gIndex_i = gRadarIndices[i];
                        int gIndex_j = gRadarIndices[j];
                        int lIndex_i = lRadarIndices[i];
                        int lIndex_j = lRadarIndices[j];
                        if(i == j) {
                            radarR(i, i) = 1;
                        }
                        else {
                            // Equation 5
                            float dist = Util::getDistance(gLocations[gIndex_i].lat(), gLocations[gIndex_i].lon(), gLocations[gIndex_j].lat(), gLocations[gIndex_j].lon(), true);
                            float h = dist / mHLengthC;
                            float rho = (1 + h) * exp(-h);
                            radarR(i, j) = rho;
                        }
                    }
                }

                float cond = arma::rcond(radarR);
                if(cond <= 0) {
                    std::stringstream ss;
                    ss << "Condition number of " << cond << " for radar values. Using raw values";
                    gridpp::util::warning(ss.str());
                    for(int e = 0; e < nEns; e++) {
                        (*output)(y, x, e) = (*field)(y, x, e); // Util::MV;
                    }
                    continue;
                }

                mattype radarRinv(lNumRadar, lNumRadar, arma::fill::zeros);
                radarRinv = arma::inv(radarR);

                for(int i = 0; i < lS; i++) {
                    int index = lLocIndices[i];
                    Rinv(i, i) = lRhos[i] / (sigma * sigma * pci[index]);
                }
                // Overwrite where we have radar pixels
                for(int i = 0; i < lNumRadar; i++) {
                    int ii = lRadarIndices[i];
                    for(int j = 0; j < lNumRadar; j++) {
                        int jj = lRadarIndices[j];
                        Rinv(ii, jj) = sqrt(lRhos[ii] * lRhos[jj]) / (sigmac * sigmac) * radarRinv(i, j);
                    }
                }
                */
            }
            else {
                abort();
            }

            // Compute C matrix
            // k x nS * nS x nS
            mattype C(nValidEns, lS);
            C = lY.t() * Rinv;

            mattype Pinv(nValidEns, nValidEns);
            float diag = 1 / delta * (nValidEns - 1);

            Pinv = C * lY + diag * arma::eye<mattype>(nValidEns, nValidEns);
            float cond = arma::rcond(Pinv);
            if(cond <= 0) {
                std::stringstream ss;
                ss << "Condition number of " << cond << ". Using raw values";
                gridpp::util::warning(ss.str());
                for(int e = 0; e < nEns; e++) {
                    output[y][x][e] = input[y][x][e]; // Util::MV;
                }
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
                print_matrix<mattype>(lY);
                std::cout << "lObs:" << std::endl;
                print_matrix<mattype>(lObs);
                std::cout << "Yhat" << std::endl;
                print_matrix<mattype>(lYhat);
            }
            eigval = sqrt(eigval);
            mattype Wcx = eigvec * arma::diagmat(eigval) * eigvec.t();
            mattype W = arma::real(Wcx);

            if(W.n_rows == 0) {
                std::stringstream ss;
                ss << "Could not find the real part of W. Using raw values.";
                gridpp::util::warning(ss.str());
                for(int e = 0; e < nEns; e++) {
                    output[y][x][e] = input[y][x][e];
                }
                continue;
            }

            // Compute PC
            mattype PC(nValidEns, lS);
            PC = P * C;

            // Compute w
            vectype w(nValidEns);
            if(mDiagnose)
                w = PC * (arma::ones<vectype>(lS));
            else
                w = PC * (lObs - lYhat);

            // Add w to W
            for(int e = 0; e < nValidEns; e++) {
                for(int e2 = 0; e2 < nValidEns; e2 ++) {
                    W(e, e2) = W(e, e2) + w(e) ;
                }
            }

            // Compute X (perturbations about model mean)
            vectype X(nValidEns);
            float total = 0;
            int count = 0;
            for(int e = 0; e < nValidEns; e++) {
                int ei = validEns[e];
                float value = input[y][x][ei];
                if(gridpp::util::is_valid(value)) {
                    X(e) = value;
                    total += value;
                    count++;
                }
                else {
                    std::cout << "Invalid value " << y << " " << x << " " << e << std::endl;
                }
            }
            float ensMean = total / count;
            for(int e = 0; e < nValidEns; e++) {
                X(e) -= ensMean;
            }

            // Write debugging information
            if(x == mX && y == mY) {
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
                std::cout << "C * lY" << std::endl;
                print_matrix<mattype>(C * lY);
                std::cout << "PC" << std::endl;
                print_matrix<mattype>(PC);
                std::cout << "W" << std::endl;
                print_matrix<mattype>(W);
                std::cout << "w" << std::endl;
                print_matrix<mattype>(w);
                std::cout << "Y:" << std::endl;
                print_matrix<mattype>(lY);
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
                    total += X(k) * W(k, e);
                }

                float currIncrement = total;

                float raw = ensMean;

                ///////////////////////////////
                // Anti-extrapolation filter //
                ///////////////////////////////
                if(!mExtrapolate) {
                    // Don't allow a final increment that is larger than any increment
                    // at station points
                    float maxInc = arma::max(lObs - (lY[e] + lYhat));
                    float minInc = arma::min(lObs - (lY[e] + lYhat));
                    if(x == mX && y == mY) {
                        std::cout << "Increments: " << maxInc << " " << minInc << " " << currIncrement << std::endl;
                    }

                    // The increment for this member. currIncrement is the increment relative to
                    // ensemble mean
                    float memberIncrement = currIncrement - X(e);
                    // Adjust increment if it gives a member increment that is outside the range
                    // of the observation increments
                    if(x == mX && y == mY) {
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
                    if(x == mX && y == mY) {
                        std::cout << "Final increment: " << currIncrement << " " << currIncrement - X(e) << std::endl;
                    }
                }
                output[y][x][ei] = ensMean + currIncrement;
            }
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
    return output;

    // return optimal_interpolation_single_member(input, bgrid.get_lats(), bgrid.get_lons(), bgrid.get_elevs(), bgrid.get_lafs(), pobs, pci, points.get_lats(), points.get_lons(), points.get_elevs(), points.get_lafs(), min_rho, hlength, vlength, wlength, max_points, elev_gradient, epsilon, output);
}
/*
int gridpp::optimal_interpolation_single_member(const vec2& input,
        const vec2& blats,
        const vec2& blons,
        const vec2& belevs,
        const vec2& blafs,
        const vec& pobs,  // gObs
        const vec& pci,   // pci
        const vec& plats,
        const vec& plons,
        const vec& pelevs,
        const vec& plafs,
        float min_rho,
        float hlength,
        float vlength,
        float wlength,
        int max_points,
        float elev_gradient,
        float epsilon,
        vec2& output) {

    gridpp::Grid bgrid(blats, blons, belevs, blafs);
    gridpp::Points points(plats, plons, pelevs, plafs);
    return optimal_interpolation_single_member(input, bgrid, pobs, pci, points, min_rho, hlength, vlength, wlength, max_points, elev_gradient, epsilon, output);
}
*/
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
