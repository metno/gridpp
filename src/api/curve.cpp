#include "gridpp.h"
#include <iostream>

using namespace gridpp;

float gridpp::apply_curve(float input, const vec& curve_ref, const vec& curve_fcst, gridpp::Extrapolation policy_below, gridpp::Extrapolation policy_above) {
    if(curve_ref.size() != curve_fcst.size())
        throw std::invalid_argument("curve_ref and curve_fcst must be the same size");
    if(curve_ref.size() == 0 || curve_ref.size() == 0)
        throw std::invalid_argument("curve_ref and curve_fcst cannot have size 0");

    int C = curve_fcst.size();
    float smallestObs  = curve_ref[0];
    float smallestFcst = curve_fcst[0];
    float largestObs   = curve_ref[C-1];
    float largestFcst  = curve_fcst[C-1];
    float output = gridpp::MV;

    // Linear interpolation within curve
    // Use equality here, because interpolate handles the case where the lowest and highest values
    // are identical. In this case, we want the average.
    if(input >= smallestFcst && input <= largestFcst) {
        output = gridpp::interpolate(input, curve_fcst, curve_ref);
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

        if (policy == gridpp::Unchanged) {
            output = input;
        }
        else {
            float slope = 1;
            if(policy == gridpp::Zero) {
                slope = 0;
            }
            else if(policy == gridpp::OneToOne || C <= 1) {
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
                    dObs  = curve_ref[C-1] - curve_ref[C-2];
                    dFcst = curve_fcst[C-1] - curve_fcst[C-2];
                }
                slope = dObs / dFcst;
            }
            else {
                throw std::invalid_argument("Unknown extrapolation policy");
            }
            output = nearestObs + slope * (input - nearestFcst);
        }
    }
    return output;
}
// QQ on a vector
vec gridpp::apply_curve(const vec& fcst, const vec& curve_ref, const vec& curve_fcst, gridpp::Extrapolation policy_below, gridpp::Extrapolation policy_above) {
    if(curve_ref.size() != curve_fcst.size())
        throw std::invalid_argument("curve_ref and curve_fcst must be the same size");
    if(curve_ref.size() == 0 || curve_ref.size() == 0)
        throw std::invalid_argument("curve_ref and curve_fcst cannot have size 0");

    int N = fcst.size();
    vec output(N, gridpp::MV);

    #pragma omp parallel for
    for(int i = 0; i < N; i++) {
        float input = fcst[i];
        output[i] = apply_curve(input, curve_ref, curve_fcst, policy_below, policy_above);
    }
    return output;
}
// QQ on a grid
vec2 gridpp::apply_curve(const vec2& fcst, const vec& curve_ref, const vec& curve_fcst, gridpp::Extrapolation policy_below, gridpp::Extrapolation policy_above) {
    int nY = fcst.size();
    if(curve_ref.size() != curve_fcst.size())
        throw std::invalid_argument("curve_ref and curve_fcst must be the same size");
    if(curve_ref.size() == 0 || curve_ref.size() == 0)
        throw std::invalid_argument("curve_ref and curve_fcst cannot have size 0");

    vec2 output(nY);
    for(int y = 0; y < nY; y++) {
        output[y] = apply_curve(fcst[y], curve_ref, curve_fcst, policy_below, policy_above);
    }
    return output;
}
// Spatially varying QQ map
vec2 gridpp::apply_curve(const vec2& fcst, const vec3& curve_ref, const vec3& curve_fcst, gridpp::Extrapolation policy_below, gridpp::Extrapolation policy_above) {
    if(!gridpp::compatible_size(curve_ref, curve_fcst))
        throw std::invalid_argument("curve_ref and curve_fcst dimension sizes mismatch");
    if(!gridpp::compatible_size(fcst, curve_ref))
        throw std::invalid_argument("Fcst and curve_ref dimension sizes mismatch");
    if(!gridpp::compatible_size(fcst, curve_fcst))
        throw std::invalid_argument("Fcst and curve_fcst dimension sizes mismatch");

    int nY = fcst.size();
    int nX = fcst[0].size();
    vec2 output(nY);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nX, gridpp::MV);
    }

    #pragma omp parallel for collapse(2)
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            float input = fcst[y][x];
            output[y][x] = apply_curve(input, curve_ref[y][x], curve_fcst[y][x], policy_below, policy_above);
        }
    }
    return output;
}
vec gridpp::monotonize_curve(vec curve_ref, vec curve_fcst, vec& output_fcst) {
    if(curve_ref.size() != curve_fcst.size())
        throw std::invalid_argument("curve_ref and curve_fcst must be the same size");
    if(curve_ref.size() == 0 || curve_ref.size() == 0)
        throw std::invalid_argument("curve_ref and curve_fcst cannot have size 0");

    bool debug = false;

    // Remove missing
    int N = curve_ref.size();
    ivec keep;
    keep.reserve(N); // Array of indices to keep (both x and y are valid)
    for(int i = 0; i < N; i++) {
        if(gridpp::is_valid(curve_fcst[i]) && gridpp::is_valid(curve_ref[i]))
            keep.push_back(i);
    }
    if(keep.size() != N) {
        vec curve_ref_copy = curve_ref;
        vec curve_fcst_copy = curve_fcst;
        curve_ref.clear();
        curve_fcst.clear();
        curve_ref.resize(keep.size(), 0);
        curve_fcst.resize(keep.size(), 0);
        for(int i = 0; i < keep.size(); i++) {
            curve_ref[i] = curve_ref_copy[keep[i]];
            curve_fcst[i] = curve_fcst_copy[keep[i]];
        }
    }
    N = curve_ref.size();
    ivec new_indices;
    new_indices.reserve(N);

    int start = 1;
    float prev = curve_fcst[0];

    /* A deviation is a part of the curve is not monotonic. This occurs when a set of x values are
       not in increasing order. Remove all points that are inside the non-increasing setion.
    */
    bool deviation = false;
    int deviation_index = 0;
    float x_min = curve_fcst[0]; // The lowest x value in the deviation section
    float x_max = curve_fcst[0]; // The highest x value in the deivation section

    new_indices.push_back(0);
    float tol = 0.1;
    for(int i = start; i < N; i++) {
        float x = curve_fcst[i];
        float y = curve_ref[i];
        assert(gridpp::is_valid(x));
        assert(gridpp::is_valid(y));
        if(deviation) {
            // Currently inside a deviation section
            if(x < x_min) {
                x_min = x;
            }
            if(x > x_max + tol) {
                // We are past the deviation section
                x_max = x;
                if(debug)
                    std::cout << "End of loop at " << i << " x=" << x << "(" << x_min << "," << x_max << ")" << std::endl;
                // Remove points inside the deviation section, that is all points above x_min
                for(int j = new_indices.size() - 1; j >= 0; j--) {
                    int index = new_indices[j];
                    if(curve_fcst[index] < x_min - tol)
                        break;
                    else {
                        if(debug)
                            std::cout << "Removing index=" << index << " x=" << curve_fcst[index] << std::endl;
                        new_indices.pop_back();
                    }
                }
                new_indices.push_back(i);
                deviation = false;
                prev = x;
                x_max = x;
            }
        }
        else {
            if(x <= prev + tol) {
                // We have a new deviation
                deviation = true;
                deviation_index = i;
                x_min = x;
                if(debug)
                    std::cout << "Detecting loop from " << i << " x=" << x << " " << x_max << std::endl;
            }
            else {
                new_indices.push_back(i);
                prev = x;
                x_max = x;
            }
        }
    }
    if(deviation) {
        if(debug)
            std::cout << "Finished in a deviation " << deviation_index << std::endl;
        // Find latest x point
        for(int j = new_indices.size() - 1; j >= 0; j--) {
            if(debug)
                std::cout << new_indices[j] << " " << deviation_index << std::endl;
            int index = new_indices[j];
            if(curve_fcst[index] >= x_min)
                new_indices.pop_back();
        }
    }

    vec output_ref(2);
    output_ref.resize(new_indices.size());
    output_fcst.clear();
    output_fcst.resize(new_indices.size());
    for(int i = 0; i < new_indices.size(); i++) {
        assert(curve_fcst.size() > new_indices[i]);
        output_ref[i] = curve_ref[new_indices[i]];
        output_fcst[i] = curve_fcst[new_indices[i]];
    }
    return output_ref;
}
