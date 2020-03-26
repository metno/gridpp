#include "gridpp.h"
#include <math.h>
vec2 gridpp::nearest(const Grid& igrid, const Grid& ogrid, vec2 ivalues) {
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec2 output(nLat);
    for(int i = 0; i < nLat; i++)
        output[i].resize(nLon);

    #pragma omp parallel for
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec indices = igrid.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            int I = indices[0];
            int J = indices[1];
            output[i][j] = ivalues[I][J];
        }
    }
    return output;
}
vec gridpp::nearest(const Grid& igrid, const Points& opoints, vec2 ivalues) {
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nP = iOutputLats.size();

    vec output(nP);

    #pragma omp parallel for
    for(int i = 0; i < nP; i++) {
        ivec indices = igrid.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        int I = indices[0];
        int J = indices[1];
        output[i] = ivalues[I][J];
    }
    return output;
}
