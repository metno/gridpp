#include "gridpp.h"

vec gridpp::apply_curve(const vec& fcst, const vec2& curve, gridpp::Extrapolation policy_below, gridpp::Extrapolation policy_above) {
    int N = fcst.size();
    int C = fcst.size();
    vec output(N, gridpp::MV);
    const vec& curve_fcst = curve[0];
    const vec& curve_ref = curve[1];
    float smallestObs  = curve_ref[0];
    float smallestFcst = curve_fcst[0];
    float largestObs   = curve_ref[C-1];
    float largestFcst  = curve_fcst[C-1];
    for(int i = 0; i < N; i++) {
        float input = fcst[i];
        // Linear interpolation within curve
        if(input > smallestFcst && input < largestFcst) {
            output[i] = gridpp::util::interpolate(input, curve_fcst, curve_ref);
        }
        // Extrapolate outside curve
        else {
            float nearestObs;
            float nearestFcst;
            gridpp::Extrapolation policy;
            if(input <= smallestFcst) {
                nearestObs  = smallestObs;
                nearestFcst = smallestFcst;
                policy = policy_below;
            }
            else {
                nearestObs  = largestObs;
                nearestFcst = largestFcst;
                policy = policy_above;
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
                    dObs  = curve_ref[1] - curve_ref[0];
                    dFcst = curve_fcst[1] - curve_fcst[0];
                }
                else {
                    dObs  = curve_ref[N-1] - curve_ref[N-2];
                    dFcst = curve_fcst[N-1] - curve_fcst[N-2];
                }
                slope = dObs / dFcst;
            }
            output[i] = nearestObs + slope * (input - nearestFcst);
        }
    }
    return output;
}
vec2 gridpp::apply_curve(const vec2& fcst, const vec2& curve, gridpp::Extrapolation policy_below, gridpp::Extrapolation policy_above) {
    int nY = fcst.size();
    int nX = fcst[0].size();
    vec2 output(nY);
    for(int y = 0; y < nY; y++) {
        output[y] = apply_curve(fcst[y], curve, policy_below, policy_above);
    }
    return output;
}
vec2 gridpp::monotonize_curve(const vec2& curve) {
    ivec new_indices;
    assert(curve.size() == 2);
    int N = curve[0].size();
    int Nmiddle = N / 2;
    new_indices.reserve(N);
    // Upper half
    int start = 0;
    float last = curve[0][start];
    bool deviation = false;
    int deviation_index = 0;
    float x_min = 0;
    float x_max = 0;
    float x_min_index = 0;
    for(int i = start; i < N; i++) {
        float x = curve[0][i];
        float y = curve[1][i];
        if(!gridpp::util::is_valid(x) || !gridpp::util::is_valid(y))
            continue;
        if(deviation) {
            if(x < x_min) {
                x_min = x;
                x_min_index = i;
            }
            if(x > x_max) {
                // Stopping criterion
                x_max = x;
                // std::cout << "End of loop at " << i << " x=" << x << "(" << x_min << "," << x_max << ")" << std::endl;
                // Find latest x point
                for(int j = new_indices.size() - 1; j >= 0; j--) {
                    int index = new_indices[j];
                    if(curve[0][index] < x_min)
                        break;
                    else {
                        // std::cout << "Removing index=" << index << " x=" << curve[0][index] << std::endl;
                        new_indices.pop_back();
                    }
                }
                new_indices.push_back(i);
            }
            deviation = false;
        }
        else {
            if(x <= last) {
                // We have a deviation
                deviation = true;
                deviation_index = i;
                x_min = x;
                x_max = x;
                x_min_index = i;
                // std::cout << "Detecting loop from " << i << " x=" << x << std::endl;
            }
            else {
                new_indices.push_back(i);
                last = x;
            }
        }
    }

    vec2 new_curve(2);
    new_curve[0].resize(new_indices.size());
    new_curve[1].resize(new_indices.size());
    for(int i = 0; i < new_indices.size(); i++) {
        assert(curve[0].size() > new_indices[i]);
        new_curve[0][i] = curve[0][new_indices[i]];
        new_curve[1][i] = curve[1][new_indices[i]];
    }
    return new_curve;
}
