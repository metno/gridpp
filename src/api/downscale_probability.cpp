#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

vec2 gridpp::downscale_probability(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec2& threshold, const ComparisonOperator& comparisonoperator) {
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec2 output(nLat);
    for(int i = 0; i < nLat; i++)
        output[i].resize(nLon);

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec indices = igrid.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            int I = indices[0];
            int J = indices[1];
            /* TODO: Implement the actual method */
        }
    }
    return output;
}