#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec gridpp::count(const Grid& grid, const Points& points, float radius) {
    int size = points.size();
    vec output(size);
    vec lats = points.get_lats();
    vec lons = points.get_lons();
    vec2 ilats = grid.get_lats();
    vec2 ilons = grid.get_lons();
    #pragma omp parallel for
    for(int i = 0; i < size; i++) {
        int num = grid.get_num_neighbours(lats[i], lons[i], radius);
        output[i] = num;
    }
    return output;
}

vec2 gridpp::count(const Grid& igrid, const Grid& ogrid, float radius) {
    ivec size = ogrid.size();
    vec2 output(size[0]);
    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    for(int i = 0; i < size[0]; i++) {
        output[i].resize(size[1], 0);
    }
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < size[0]; i++) {
        for(int j = 0; j < size[1]; j++) {
            int num = igrid.get_num_neighbours(ilats[i][j], ilons[i][j], radius);
            output[i][j] = num;
        }
    }
    return output;
}

vec2 gridpp::count(const Points& points, const Grid& grid, float radius) {
    ivec size = grid.size();
    vec2 output(size[0]);
    vec lats = points.get_lats();
    vec lons = points.get_lons();
    vec2 ilats = grid.get_lats();
    vec2 ilons = grid.get_lons();
    for(int i = 0; i < size[0]; i++) {
        output[i].resize(size[1], 0);
    }
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < size[0]; i++) {
        for(int j = 0; j < size[1]; j++) {
            int num = points.get_num_neighbours(ilats[i][j], ilons[i][j], radius);
            output[i][j] = num;
        }
    }
    return output;
}

vec gridpp::count(const Points& ipoints, const Points& opoints, float radius) {
    int size = opoints.size();
    vec output(size);
    vec ilats = ipoints.get_lats();
    vec ilons = ipoints.get_lons();
    vec olats = opoints.get_lats();
    vec olons = opoints.get_lons();
    #pragma omp parallel for
    for(int i = 0; i < size; i++) {
        int num = ipoints.get_num_neighbours(ilats[i], ilons[i], radius);
        output[i] = num;
    }
    return output;
}
