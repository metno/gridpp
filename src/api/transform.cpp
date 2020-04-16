#include "gridpp.h"
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
float gridpp::Log::forward(float value) const {
    return log(value);
}
float gridpp::Log::backward(float value) const {
    return exp(value);
}
gridpp::BoxCox::BoxCox(float threshold) : mThreshold(threshold) {

}
float gridpp::BoxCox::forward(float value) const {
   if(value <= 0)
      value = 0;
   if(mThreshold == 0)
      return log(value);
   else
      return (pow(value, mThreshold) - 1) / mThreshold;
}
float gridpp::BoxCox::backward(float value) const {
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
float gridpp::Identity::forward(float value) const {
    return value;
}
float gridpp::Identity::backward(float value) const {
    return value;
}
