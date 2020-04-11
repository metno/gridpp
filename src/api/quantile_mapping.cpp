#include "gridpp.h"

namespace {
    void separate(const vec& parameters, vec& ref, vec& fcst);
}

vec2 gridpp::quantile_mapping(const gridpp::Grid& grid, const vec2& input, gridpp::Extrapolation policy, const gridpp::Parameters& parameters) {
    int Y = input.size();
    int X = input[0].size();
    vec2 lats = grid.get_lats();
    vec2 lons = grid.get_lons();
    vec2 altitude = grid.get_elevs();
    vec2 land_area_fraction = grid.get_lafs();
    vec2 output(Y);
    for(int y = 0; y < Y; y++) {
        output[y].resize(X, gridpp::MV);
        for(int x = 0; x < X; x++) {
            vec par = parameters.get(lats[y][x], lons[y][x], altitude[y][x], land_area_fraction[y][x]);
            vec fcst, ref;
            ::separate(par, ref, fcst);
            output[y][x] = gridpp::quantile_mapping(input[y][x], ref, fcst, policy);
        }
    }
    return output;
}
vec2 gridpp::quantile_mapping(const vec2& input, const vec& x, const vec& y, gridpp::Extrapolation policy) {
    gridpp::util::not_implemented_error();
    return vec2();
}
vec gridpp::quantile_mapping(const vec& input, const vec& x, const vec& y, gridpp::Extrapolation policy) {
    gridpp::util::not_implemented_error();
    return vec();
}
float gridpp::quantile_mapping(float input, const vec& ref, const vec& fcst, gridpp::Extrapolation policy) {
    float output = gridpp::MV;
    if(gridpp::util::is_valid(input)) {
        float smallestObs  = ref[0];;
        float smallestFcst = fcst[0];
        int N = ref.size();
        float largestObs   = ref[N-1];
        float largestFcst  = fcst[N-1];

        // Linear interpolation within curve
        if(input > smallestFcst && input < largestFcst) {
            output = gridpp::util::interpolate(input, fcst, ref);
        }
        // Extrapolate outside curve
        else {
            float nearestObs;
            float nearestFcst;
            if(input <= smallestFcst) {
                nearestObs  = smallestObs;
                nearestFcst = smallestFcst;
            }
            else {
                nearestObs  = largestObs;
                nearestFcst = largestFcst;
            }
            float slope = 1;
            if(policy == gridpp::ExtrapolationZero) {
                slope = 0;
            }
            if(policy == gridpp::ExtrapolationOneToOne || N <= 1) {
                slope = 1;
            }
            else if(policy == gridpp::ExtrapolationMeanSlope) {
                float dObs  = largestObs - smallestObs;
                float dFcst = largestFcst - smallestFcst;
                slope = dObs / dFcst;
            }
            else if(policy == gridpp::ExtrapolationNearestSlope) {
                float dObs;
                float dFcst;
                if(input <= smallestFcst) {
                    dObs  = ref[1] - ref[0];
                    dFcst = fcst[1] - fcst[0];
                }
                else {
                    dObs  = ref[N-1] - ref[N-2];
                    dFcst = fcst[N-1] - fcst[N-2];
                }
                slope = dObs / dFcst;
            }
            output = nearestObs + slope * (input - nearestFcst);
        }
    }
    return output;
}

namespace {
    void separate(const vec& parameters, vec& ref, vec& fcst) {
       ref.clear();
       fcst.clear();
       int N = parameters.size() / 2;
       ref.resize(N, gridpp::MV);
       fcst.resize(N, gridpp::MV);
       for(int i = 0; i < N; i++) {
          ref[i] = parameters[2*i];
          fcst[i] = parameters[2*i+1];
       }
    }
}
