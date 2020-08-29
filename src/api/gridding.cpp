#include "gridpp.h"
#include <iostream>

vec2 gridpp::gridding(const Grid& grid, const Points& points, const vec& values, float radius, int min_num, gridpp::Statistic statistic) {
    int Y = grid.size()[0];
    int X = grid.size()[1];
    vec2 lats = grid.get_lats();
    vec2 lons = grid.get_lons();
    vec2 output = gridpp::util::init_vec2(Y, X);

    // Compute the statistic of the poin values at each gridpoint
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            ivec I = points.get_neighbours(lats[y][x], lons[y][x], radius);
            if(min_num > 0 && I.size() >= min_num) {
                vec curr(I.size());
                for(int i = 0; i < I.size(); i++) {
                    curr[i] = values[I[i]];
                }
                output[y][x] = gridpp::util::calc_statistic(curr, statistic);
                if(output[y][x] >= 1000) {
                    std::cout << y << " " << x << " " << output[y][x] << std::endl;
                    for(int i = 0; i < I.size(); i++) {
                        std::cout << "   " << curr[i] << std::endl;
                    }
                }
            }
        }
    }
    return output;
}
