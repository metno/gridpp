#include "gridpp.h"
#include <math.h>
#include <algorithm>
#include <armadillo>
#include <assert.h>

int gridpp::optimal_interpolation_single_member(const vec2& input,
        const gridpp::KDTree& btree,
        const vec& pobs,  // gObs
        const vec& pci,   // gCi
        const gridpp::KDTree& ptree,
        float minRho,
        float hlength,
        float vlength,
        float wmin,
        float maxElevDiff,
        bool landOnly,
        int maxLocations,
        float elevGradient,
        float epsilon,
        vec2& output) {
    return optimal_interpolation_single_member(input, btree.get_lats_2d(), btree.get_lons_2d(), btree.get_elevs_2d(), btree.get_lafs_2d(), pobs, pci, ptree.get_lats(), ptree.get_lons(), ptree.get_elevs(), ptree.get_lafs(), minRho, hlength, vlength, wmin, maxElevDiff, landOnly, maxLocations, elevGradient, epsilon, output);
}
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
        float maxElevDiff,
        bool landOnly,
        int maxLocations,
        float elevGradient,
        float epsilon,
        vec2& output) {

    std::cout << "Starting" << std::endl;

    int nY = input.size();
    int nX = input[0].size();
    int nS = pobs.size();

    // Check input data
    check_vec(blats, nY, nX);
    check_vec(blons, nY, nX);
    check_vec(belevs, nY, nX);
    check_vec(blafs, nY, nX);
    check_vec(pci, nS);
    check_vec(plats, nS);
    check_vec(plons, nS);
    check_vec(pelevs, nS);
    check_vec(plafs, nS);

    // Prepare output matrix
    output.resize(nY);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nX);
    }

    float MV = -999;

    // Estimate the grid spacing
    float gridSize = getDistance(blats[0][0], blons[0][0], blats[1][0], blons[1][0], false);
    std::stringstream ss;
    ss << "Estimated grid size: " << gridSize << " m" << std::endl;
    ss << "Number of observations: " << nS << std::endl;
    ss << "Number of gridpoints: " << nY << " " << nX;
    debug(ss.str());

    // Loop over each observation, find the nearest gridpoint and place the obs into all gridpoints
    // in the vicinity of the nearest neighbour. This is only meant to be an approximation, but saves
    // considerable time instead of doing a loop over each grid point and each observation.

    // Store the indicies (into the gLocations array) that a gridpoint has available
    std::vector<std::vector<std::vector<int> > > gLocIndices; // Y, X, obs indices
    std::vector<float> pYi(nS, MV);
    std::vector<float> pXi(nS, MV);
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

    std::cout << "Computing search tree..." << std::endl;
    gridpp::KDTree searchTree(blats0, blons0);
    std::cout << "Done computing search tree" << std::endl;

    // For each gridpoint, find which observations are relevant. Parse the observations and only keep
    // those that pass certain checks
    for(int i = 0; i < nS; i++) {
        int index = searchTree.get_nearest_neighbour(plats[i], plons[i]);
        // TODO: Check
        int X = index / nY;
        int Y = index % nY;
        pYi[i] = Y;
        pXi[i] = X;

        // Check if the elevation of the station roughly matches the reference grid elevation
        bool hasValidElev = !isValid(maxElevDiff) || isValid(pelevs[i]);
        if(isValid(maxElevDiff) && hasValidElev) {
            float elevDiff = abs(pelevs[i] - belevs[Y][X]);
            hasValidElev = elevDiff < maxElevDiff;
        }
        if(hasValidElev) {
            // Don't include an observation if it is in the ocean and landOnly=1
            bool wrongLaf = isValid(plafs[i]) && landOnly && plafs[i] == 0;
            if(!wrongLaf) {
                for(int y = std::max(0, Y - gridpointRadius); y < std::min(nY, Y + gridpointRadius); y++) {
                    for(int x = std::max(0, X - gridpointRadius); x < std::min(nX, X + gridpointRadius); x++) {
                        gLocIndices[y][x].push_back(i);
                    }
                }
            }
        }
    }
    std::cout << "Done assigning observations to gridpoints" << std::endl;

    // Transform the background
    // #pragma omp parallel for
    // for(int x = 0; x < nX; x++) {
    //     for(int y = 0; y < nY; y++) {
    //         float value = input[y][x];
    //         if(isValid(value))
    //             input[y][x] = gridpp::transform(value);
    //     }
    // }

    // Compute Y
    vec gY(nS);
    for(int i = 0; i < nS; i++) {
        float elevCorr = 0;
        if(isValid(elevGradient) && elevGradient != 0) {
            float nnElev = belevs[pYi[i]][pXi[i]];
            assert(isValid(nnElev));
            assert(isValid(pelevs[i]));
            float elevDiff = pelevs[i] - nnElev;
            elevCorr = elevGradient * elevDiff;
        }
        float value = input[pYi[i]][pXi[i]];
        // if(isValid(value)) {
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
            std::vector<int> lLocIndices;
            lLocIndices.reserve(lLocIndices0.size());
            std::vector<std::pair<float,int> > lRhos0;
            lRhos0.reserve(lLocIndices0.size());
            for(int i = 0; i < lLocIndices0.size(); i++) {
                int index = lLocIndices0[i];
                float hdist = getDistance(plats[index], plons[index], lat, lon, true);
                float vdist = MV;
                if(isValid(pelevs[index] && isValid(elev)))
                    vdist = pelevs[index] - elev;
                float lafdist = 0;
                if(isValid(plafs[index]) && isValid(laf))
                    lafdist = plafs[index] - laf;
                float rho = gridpp::calcRho(hdist, vdist, lafdist, hlength, vlength, wmin);
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
            if(lRhos0.size() > maxLocations) {
                // If sorting is enabled and we have too many locations, then only keep the best ones based on rho.
                // Otherwise, just use the last locations added
                lRhos = arma::vec(maxLocations);
                std::sort(lRhos0.begin(), lRhos0.end(), gridpp::sort_pair_first<float,int>());
                for(int i = 0; i < maxLocations; i++) {
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
                float hdist = getDistance(plats[index], plons[index], lat, lon, true);
                float vdist = MV;
                if(isValid(pelevs[index] && isValid(elev)))
                    vdist = pelevs[index] - elev;
                float lafdist = 0;
                if(isValid(plafs[index]) && isValid(laf))
                    lafdist = plafs[index] - laf;
                float rho = calcRho(hdist, vdist, lafdist, hlength, vlength, wmin);
                lG(0, i) = rho;
                for(int j = 0; j < lS; j++) {
                    int index_j = lLocIndices[j];
                    float hdist = getDistance(plats[index], plons[index], plats[index_j], plons[index_j], true);
                    float vdist = MV;
                    if(isValid(pelevs[index] && isValid(pelevs[index_j])))
                        vdist = pelevs[index] - pelevs[index_j];
                    float lafdist = 0;
                    if(isValid(plafs[index]) && isValid(laf))
                        lafdist = plafs[index] - plafs[index_j];

                    lP(i, j) = gridpp::calcRho(hdist, vdist, lafdist, hlength, vlength, wmin);
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
            if(isValid(value))
                (*output)(y, x, e) = invTransform(value);
        }
    }
    */
    return 0;
}
float gridpp::calcRho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wmin) {
   float h = (iHDist/hlength);
   float rho = exp(-0.5 * h * h);
   if(isValid(vlength)) {
      if(!isValid(iVDist)) {
         rho = 0;
      }
      else {
         float v = (iVDist/vlength);
         rho *= exp(-0.5 * v * v);
      }
   }
   if(isValid(wmin)) {
      rho *= 1 - (1 - wmin) * std::abs(iLDist);
   }
   return rho;
}
float gridpp::getDistance(float lat1, float lon1, float lat2, float lon2, bool approx) {
   if(!isValid(lat1) || !isValid(lat2) ||
      !isValid(lon1) || !isValid(lon2)) {
      return -999;
   }
   if(!(fabs(lat1) <= 90 && fabs(lat2) <= 90 && fabs(lon1) <= 360 && fabs(lon2) <= 360)) {
      std::stringstream ss;
      ss  <<" Cannot calculate distance, invalid lat/lon: (" << lat1 << "," << lon1 << ") (" << lat2 << "," << lon2 << ")";
      error(ss.str());
   }

   if(lat1 == lat2 && lon1 == lon2)
      return 0;

   double lat1r = deg2rad(lat1);
   double lat2r = deg2rad(lat2);
   double lon1r = deg2rad(lon1);
   double lon2r = deg2rad(lon2);
   // Calculate distance according to: http://www.movable-type.co.uk/scripts/latlong.html
   double radiusEarth = 6.378137e6;
   if(approx) {
      float dx2 = pow(cos((lat1r+lat2r)/2),2)*(lon1r-lon2r)*(lon1r-lon2r);
      float dy2 = (lat1r-lat2r)*(lat1r-lat2r);
      return radiusEarth*sqrt(dx2+dy2);
   }
   else {
      double ratio = cos(lat1r)*cos(lon1r)*cos(lat2r)*cos(lon2r)
                   + cos(lat1r)*sin(lon1r)*cos(lat2r)*sin(lon2r)
                   + sin(lat1r)*sin(lat2r);
      double dist = acos(ratio)*radiusEarth;
      return (float) dist;
   }
}
bool gridpp::isValid(float iValue) {
    return !std::isnan(iValue) && !std::isinf(iValue) && iValue != -999;
}
float gridpp::deg2rad(float deg) {
   float pi = 3.14159265;
   return (deg * pi / 180);
}
float gridpp::rad2deg(float rad) {
   float pi = 3.14159265;
   return (rad * 180 / pi);
}
// Set up convenient functions for debugging in gdb
// template<class Matrix>
// void print_matrix(Matrix matrix) {
//     matrix.print(std::cout);
// }

// template void print_matrix<CalibratorOi::mattype>(CalibratorOi::mattype matrix);
// template void print_matrix<CalibratorOi::cxtype>(CalibratorOi::cxtype matrix);
void gridpp::debug(std::string string) {
    std::cout << string << std::endl;
}

void gridpp::error(std::string string) {
    std::cout << string << std::endl;
}
void gridpp::check_vec(vec2 input, int Y, int X) {
    assert(input.size() == Y);
    for(int i = 0; i < input.size(); i++) {
        assert(input[i].size() == X);
        for(int j = 0; j < input[i].size(); j++) {
            assert(isValid(input[i][j]));
        }
    }
}
void gridpp::check_vec(vec input, int S) {
    assert(input.size() == S);
    for(int i = 0; i < input.size(); i++) {
        assert(isValid(input[i]));
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
