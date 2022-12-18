#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec gridpp::count(const Grid& grid, const Points& points, float radius) {
    int size = points.size();
    vec output(size);
    vec olats = points.get_lats();
    vec olons = points.get_lons();
    #pragma omp parallel for
    for(int i = 0; i < size; i++) {
        int num = grid.get_num_neighbours(olats[i], olons[i], radius);
        output[i] = num;
    }
    return output;
}

vec2 gridpp::count(const Grid& igrid, const Grid& ogrid, float radius) {
    ivec size = ogrid.size();
    vec2 output(size[0]);
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    for(int i = 0; i < size[0]; i++) {
        output[i].resize(size[1], 0);
    }
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < size[0]; i++) {
        for(int j = 0; j < size[1]; j++) {
            int num = igrid.get_num_neighbours(olats[i][j], olons[i][j], radius);
            output[i][j] = num;
        }
    }
    return output;
}

vec2 gridpp::count(const Points& points, const Grid& grid, float radius) {
    ivec size = grid.size();
    vec2 output(size[0]);
    vec2 olats = grid.get_lats();
    vec2 olons = grid.get_lons();
    for(int i = 0; i < size[0]; i++) {
        output[i].resize(size[1], 0);
    }
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < size[0]; i++) {
        for(int j = 0; j < size[1]; j++) {
            int num = points.get_num_neighbours(olats[i][j], olons[i][j], radius);
            output[i][j] = num;
        }
    }
    return output;
}

vec gridpp::count(const Points& ipoints, const Points& opoints, float radius) {
    int size = opoints.size();
    vec output(size);
    vec olats = opoints.get_lats();
    vec olons = opoints.get_lons();
    #pragma omp parallel for
    for(int i = 0; i < size; i++) {
        int num = ipoints.get_num_neighbours(olats[i], olons[i], radius);
        output[i] = num;
    }
    return output;
}
