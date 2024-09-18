#include <iostream>
#include "gridpp.h"

using namespace gridpp;

// These two are included because of SWIG. See note about SWIG requirement in gridpp.h
float gridpp::Transform::forward(float value) const {
    return -1;
}
float gridpp::Transform::backward(float value) const {
    return -1;
}
vec gridpp::Transform::forward(const vec& input) const {
    int Y = input.size();
    vec output(Y, gridpp::MV);
    for(int y = 0; y < Y; y++) {
        output[y] = forward(input[y]);
    }
    return output;
}
vec2 gridpp::Transform::forward(const vec2& input) const {
    int Y = input.size();
    vec2 output(Y);
    for(int y = 0; y < Y; y++) {
        int X = input[y].size();
        output[y].resize(X, gridpp::MV);
        for(int x = 0; x < X; x++) {
            output[y][x] = forward(input[y][x]);
        }
    }
    return output;
}
vec3 gridpp::Transform::forward(const vec3& input) const {
    int Y = input.size();
    vec3 output(Y);
    for(int y = 0; y < Y; y++) {
        int X = input[y].size();
        output[y].resize(X);
        for(int x = 0; x < X; x++) {
            int E = input[y][x].size();
            output[y][x].resize(E, gridpp::MV);
            for(int e = 0; e < E; e++) {
                output[y][x][e] = forward(input[y][x][e]);
            }
        }
    }
    return output;
}
vec gridpp::Transform::backward(const vec& input) const {
    int Y = input.size();
    vec output(Y, gridpp::MV);
    for(int y = 0; y < Y; y++) {
        output[y] = backward(input[y]);
    }
    return output;
}
vec2 gridpp::Transform::backward(const vec2& input) const {
    int Y = input.size();
    vec2 output(Y);
    for(int y = 0; y < Y; y++) {
        int X = input[y].size();
        output[y].resize(X, gridpp::MV);
        for(int x = 0; x < X; x++) {
            output[y][x] = backward(input[y][x]);
        }
    }
    return output;
}
vec3 gridpp::Transform::backward(const vec3& input) const {
    int Y = input.size();
    vec3 output(Y);
    for(int y = 0; y < Y; y++) {
        int X = input[y].size();
        output[y].resize(X);
        for(int x = 0; x < X; x++) {
            int E = input[y][x].size();
            output[y][x].resize(E, gridpp::MV);
            for(int e = 0; e < E; e++) {
                output[y][x][e] = backward(input[y][x][e]);
            }
        }
    }
    return output;
}
float gridpp::Log::forward(float value) const {
    if(gridpp::is_valid(value))
        return log(value);
    else
        return gridpp::MV;
}
float gridpp::Log::backward(float value) const {
    if(gridpp::is_valid(value))
        return exp(value);
    else
        return gridpp::MV;
}
gridpp::BoxCox::BoxCox(float threshold) : mThreshold(threshold) {

}
float gridpp::BoxCox::forward(float value) const {
    if(!gridpp::is_valid(value))
        return gridpp::MV;
    if(value <= 0)
        value = 0;
    if(mThreshold == 0)
        return log(value);
    else
        return (pow(value, mThreshold) - 1) / mThreshold;
}
float gridpp::BoxCox::backward(float value) const {
    if(!gridpp::is_valid(value))
        return gridpp::MV;
   float rValue = 0;
   if(mThreshold == 0)
      rValue = exp(value);
   else {
      if(value < -1.0 / mThreshold) {
         value = -1.0 / mThreshold;
      }
      rValue = pow(1 + mThreshold * value, 1 / mThreshold);
   }
   if(rValue <= 0)
      rValue = 0;
   return rValue;
}
float gridpp::StartedBoxCox::forward(float value) const {
    if(!gridpp::is_valid(value) || mThreshold <= 0)
        return gridpp::MV;
    if(value <= 0)
        value = 0;
    if(value <= 1)
        return value;
    else
        return (1 + (((pow(value, mThreshold)) - 1) / mThreshold));
}
float gridpp::StartedBoxCox::backward(float value) const {
    if(!gridpp::is_valid(value) || mThreshold <= 0)
        return gridpp::MV;
    float rValue = 0;
    if(value <= 1)
        rValue = value;
    else 
        rValue = pow( 1 + mThreshold * (value-1), 1 / mThreshold);
    if(rValue <= 0)
       rValue = 0;
    return rValue;
}
gridpp::Gamma::Gamma(float shape, float scale, float tolerance) : m_gamma_dist(1, 1), m_norm_dist(), m_tolerance(tolerance) {
    // Initialize the gamma distribution to something that works, and then overwrite it so that
    // we can check for argument errors gracefully
    if(!gridpp::is_valid(shape) || shape <= 0)
        throw std::invalid_argument("Shape parameter must be > 0 in the gamma distribution");
    if(!gridpp::is_valid(scale) || scale <= 0)
        throw std::invalid_argument("Scale parameter must be > 0 in the gamma distribution");
    if(!gridpp::is_valid(tolerance) || tolerance < 0)
        throw std::invalid_argument("Tolerance must be >= 0 in the gamma distribution");
    m_gamma_dist = boost::math::gamma_distribution<> (shape, scale);
}
float gridpp::Gamma::forward(float value) const {
    if(!gridpp::is_valid(value))
        return gridpp::MV;
    float cdf = boost::math::cdf(m_gamma_dist, value + m_tolerance);
    float result = boost::math::quantile(m_norm_dist, cdf);
    return result;
}
float gridpp::Gamma::backward(float value) const {
    if(!gridpp::is_valid(value))
        return gridpp::MV;
   float cdf = boost::math::cdf(m_norm_dist, value);
   float result = boost::math::quantile(m_gamma_dist, cdf) - m_tolerance;
   return result;
}
float gridpp::Identity::forward(float value) const {
    return value;
}
float gridpp::Identity::backward(float value) const {
    return value;
}
