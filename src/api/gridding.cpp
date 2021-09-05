#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec2 gridpp::gridding(const Grid& grid, const Points& points, const vec& values, float radius, int min_num, gridpp::Statistic statistic) {
    if(!gridpp::compatible_size(points, values))
        throw std::invalid_argument("Points size is not the same as values");
    int Y = grid.size()[0];
    int X = grid.size()[1];
    vec2 lats = grid.get_lats();
    vec2 lons = grid.get_lons();
    vec2 output = gridpp::init_vec2(Y, X);

    // Compute the statistic of the point values at each gridpoint
    #pragma omp parallel for collapse(2)
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            ivec I = points.get_neighbours(lats[y][x], lons[y][x], radius);
            if(min_num <= 0 || I.size() >= min_num) {
                vec curr(I.size());
                for(int i = 0; i < I.size(); i++) {
                    curr[i] = values[I[i]];
                }
                output[y][x] = gridpp::calc_statistic(curr, statistic);
            }
        }
    }
    return output;
}
vec2 gridpp::gridding_nearest(const Grid& grid, const Points& points, const vec& values, int min_num, gridpp::Statistic statistic) {
    if(!gridpp::compatible_size(points, values))
        throw std::invalid_argument("Points size is not the same as values");
    int Y = grid.size()[0];
    int X = grid.size()[1];
    vec lats = points.get_lats();
    vec lons = points.get_lons();
    vec3 temp(Y);
    for(int y = 0; y < Y; y++) {
        temp[y].resize(X);
    }

    int S = values.size();

    // Add value to the nearest grid point in the grid
    // Not parallelizable, because two threads might be writing to the same piece of memory
    for(int s = 0; s < S; s++) {
        ivec indices = grid.get_nearest_neighbour(lats[s], lons[s]);
        int iY = indices[0];
        int iX = indices[1];
        temp[iY][iX].push_back(values[s]);
    }

    vec2 output = gridpp::init_vec2(Y, X, gridpp::MV);
    #pragma omp parallel for collapse(2)
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            if(temp[y][x].size() > 0) {
                if(min_num <=0 || temp[y][x].size() >= min_num)
                    output[y][x] = gridpp::calc_statistic(temp[y][x], statistic);
            }
        }
    }
    return output;
}
