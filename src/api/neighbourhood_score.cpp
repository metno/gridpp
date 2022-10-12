#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec2 gridpp::neighbourhood_score(const Grid& grid, const Points& points, const vec2& fcst, const vec& ref, int half_width, gridpp::Metric metric, float threshold) {

    if(!gridpp::compatible_size(grid, fcst)) {
        throw std::invalid_argument("Grid size is not the same as forecast values");
    }

    if(half_width <= 0) {
        throw std::invalid_argument("half_width must be greater than 0");
    }

    int nY = fcst.size();
    int nX = fcst[0].size();

    vec2 a = gridpp::init_vec2(nY, nX, 0);
    vec2 b = gridpp::init_vec2(nY, nX, 0);
    vec2 c = gridpp::init_vec2(nY, nX, 0);
    vec2 d = gridpp::init_vec2(nY, nX, 0);

    // Gridding of observations on the the forecast grid
    vec2 ref_grid = gridpp::gridding_nearest(grid, points, ref, 1, gridpp::Mean);

    // Compute the 4 contingency values
    #pragma omp parallel for collapse(2)
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            if(gridpp::is_valid(ref_grid[y][x]) && gridpp::is_valid(fcst[y][x])) {
                if(fcst[y][x] > threshold) {
                    a[y][x] = ref_grid[y][x] > threshold;
                    b[y][x] = ref_grid[y][x] <= threshold;
                }
                else {
                    c[y][x] = ref_grid[y][x] > threshold;
                    d[y][x] = ref_grid[y][x] <= threshold;
                }
            }
        }
    }

    // Apply neighbourhood method on contingency values
    vec2 a_hood = gridpp::neighbourhood(a, half_width, gridpp::Mean);
    vec2 b_hood = gridpp::neighbourhood(b, half_width, gridpp::Mean);
    vec2 c_hood = gridpp::neighbourhood(c, half_width, gridpp::Mean);
    vec2 d_hood = gridpp::neighbourhood(d, half_width, gridpp::Mean);

    // Compute score in neighbourhood
    vec2 output = gridpp::init_vec2(nY, nX, gridpp::Mean);
    #pragma omp parallel for collapse(2)
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            output[y][x] = gridpp::calc_score(a_hood[y][x], b_hood[y][x], c_hood[y][x], d_hood[y][x], metric);
            // std::cout << y << " " << x << " " << a_hood[y][x] << std::endl;
        }
    }
    return output;
}
