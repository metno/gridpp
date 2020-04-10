#include "gridpp.h"
#include <cmath>

float gridpp::wind_speed(float xwind, float ywind) {
    return sqrt(xwind * xwind + ywind * ywind);
}
vec gridpp::wind_speed(const vec& xwind, const vec& ywind) {
    int N = xwind.size();
    vec values(N, gridpp::MV);
    for(int n = 0; n < N; n++) {
        values[n] = wind_speed(xwind[n], ywind[n]);
    }
    return values;
}
vec2 gridpp::wind_speed(const vec2& xwind, const vec2& ywind) {
    int Y = xwind.size();
    vec2 values(Y);
    if(Y > 0) {
        int X = xwind[0].size();
        for(int y = 0; y < Y; y++) {
            values[y] = wind_speed(xwind[y], ywind[y]);
        }
    }
    return values;
}
