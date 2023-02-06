#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

vec2 gridpp::nearest(const Grid& igrid, const Grid& ogrid, const vec2& ivalues) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec2 output = gridpp::init_vec2(nLat, nLon, gridpp::MV);
    if(igrid.size()[0] == 0 || igrid.size()[1] == 0)
        return output;

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec indices = igrid.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            assert(indices.size() > 0);
            int I = indices[0];
            int J = indices[1];
            output[i][j] = ivalues[I][J];
        }
    }
    return output;
}
vec3 gridpp::nearest(const Grid& igrid, const Grid& ogrid, const vec3& ivalues) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nTime = ivalues.size();
    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();
;
    vec3 output = gridpp::init_vec3(nTime, nLat, nLon, gridpp::MV);
    if(igrid.size()[0] == 0 || igrid.size()[1] == 0)
        return output;

    // Fetch the nearest neighbours first, because we do not want time to be the inner loop. For
    // large number of times, the memory access gets slow when time is the inner loop.
    ivec2 I = gridpp::init_ivec2(nLat, nLon, -1);
    ivec2 J = gridpp::init_ivec2(nLat, nLon, -1);

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec indices = igrid.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            I[i][j] = indices[0];
            J[i][j] = indices[1];
        }
    }

    #pragma omp parallel for collapse(2)
    for(int t = 0; t < nTime; t++) {
        for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
                int II = I[i][j];
                int JJ = J[i][j];
                output[t][i][j] = ivalues[t][II][JJ];
            }
        }
    }
    return output;
}

vec2 gridpp::nearest(const Points& ipoints, const Grid& ogrid, const vec& ivalues) {
    if(ipoints.size() != ivalues.size())
        throw std::invalid_argument("Points size is not the same as values");
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec2 output = gridpp::init_vec2(nLat, nLon, gridpp::MV);
    if(ipoints.size() == 0)
        return output;

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            int index = ipoints.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            output[i][j] = ivalues[index];
        }
    }
    return output;
}
vec3 gridpp::nearest(const Points& ipoints, const Grid& ogrid, const vec2& ivalues) {
    if(!gridpp::compatible_size(ipoints, ivalues)) {
        std::stringstream ss;
        ss << "Input points (" << ipoints.size() << ") and values (" << ivalues.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nTime = ivalues.size();
    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec3 output = gridpp::init_vec3(nTime, nLat, nLon, gridpp::MV);
    if(ipoints.size() == 0)
        return output;

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            int index = ipoints.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            for(int t = 0; t < nTime; t++) {
                output[t][i][j] = ivalues[t][index];
            }
        }
    }
    return output;
}

vec gridpp::nearest(const Grid& igrid, const Points& opoints, const vec2& ivalues) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nPoints = iOutputLats.size();

    vec output(nPoints, gridpp::MV);
    if(igrid.size()[0] == 0 || igrid.size()[1] == 0)
        return output;

    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        ivec indices = igrid.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        int I = indices[0];
        int J = indices[1];
        output[i] = ivalues[I][J];
    }
    return output;
}
vec2 gridpp::nearest(const Grid& igrid, const Points& opoints, const vec3& ivalues) {
    if(!gridpp::compatible_size(igrid, ivalues))
       throw std::invalid_argument("Grid size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nPoints = iOutputLats.size();
    int nTime = ivalues.size();

    vec2 output = gridpp::init_vec2(nTime, nPoints, gridpp::MV);
    if(igrid.size()[0] == 0 || igrid.size()[1] == 0)
        return output;

    // Fetch the nearest neighbours first, because we do not want time to be the inner loop. For
    // large number of times, the memory access gets slow when time is the inner loop.
    ivec I(nPoints, -1);
    ivec J(nPoints, -1);
    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        ivec indices = igrid.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        I[i] = indices[0];
        J[i] = indices[1];
    }

    #pragma omp parallel for collapse(2)
    for(int t = 0; t < nTime; t++) {
        for(int i = 0; i < nPoints; i++) {
            output[t][i] = ivalues[t][I[i]][J[i]];
        }
    }
    return output;
}
vec gridpp::nearest(const Points& ipoints, const Points& opoints, const vec& ivalues) {
    if(ipoints.size() != ivalues.size())
        throw std::invalid_argument("Points size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nPoints = iOutputLats.size();

    vec output(nPoints, gridpp::MV);
    if(ipoints.size() == 0)
        return output;

    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        int index = ipoints.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        assert(index != -1);
        assert(index < ivalues.size());
        output[i] = ivalues[index];
    }
    return output;
}
vec2 gridpp::nearest(const Points& ipoints, const Points& opoints, const vec2& ivalues) {
    if(!gridpp::compatible_size(ipoints, ivalues)) {
        std::stringstream ss;
        ss << "Input points (" << ipoints.size() << ") and values (" << ivalues[0].size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nPoints = iOutputLats.size();
    int nTime = ivalues.size();

    vec2 output = gridpp::init_vec2(nTime, nPoints, gridpp::MV);
    if(ipoints.size() == 0)
        return output;

    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        int index = ipoints.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        for(int t = 0; t < nTime; t++) {
            output[t][i] = ivalues[t][index];
        }
    }
    return output;
}
