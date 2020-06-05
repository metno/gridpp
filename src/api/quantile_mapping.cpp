#include "gridpp.h"

namespace {
    void separate(const vec& parameters, vec& ref, vec& fcst);
}

vec2 gridpp::quantile_mapping(const vec2& input, const vec& x, const vec& y, gridpp::Extrapolation policy) {
    gridpp::util::not_implemented_error();
    return vec2();
}
vec gridpp::quantile_mapping(const vec& input, const vec& x, const vec& y, gridpp::Extrapolation policy) {
    gridpp::util::not_implemented_error();
    return vec();
}
float gridpp::quantile_mapping(float input, const vec& fcst, const vec& ref, gridpp::Extrapolation policy) {
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
            if(policy == gridpp::Zero) {
                slope = 0;
            }
            if(policy == gridpp::OneToOne || N <= 1) {
                slope = 1;
            }
            else if(policy == gridpp::MeanSlope) {
                float dObs  = largestObs - smallestObs;
                float dFcst = largestFcst - smallestFcst;
                slope = dObs / dFcst;
            }
            else if(policy == gridpp::NearestSlope) {
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
