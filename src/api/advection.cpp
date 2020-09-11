#include "gridpp.h"
#include <iostream>


vec2 gridpp::advection_implicit(const vec2& y_dist, const vec2& x_dist, float dt, ivec2& y_coord, ivec2& x_coord) {
    int Y = y_dist.size();
    int X = y_dist[0].size();
    y_coord.clear();
    x_coord.clear();
    y_coord.resize(Y);
    x_coord.resize(Y);
    for(int y = 0; y < Y; y++) {
        y_coord[y].resize(Y);
        x_coord[y].resize(Y);
        for(int x = 0; x < X; x++) {
            y_coord[y][x] = gridpp::MV;
            x_coord[y][x] = gridpp::MV;
        }
    }
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            int l = int(y + dt * y_dist[y][x]);
            int k = int(x + dt * x_dist[y][x]);
            if(k >= 0 && k < X && l >= 0 && l < Y) {
                y_coord[l][k] = y;
                x_coord[l][k] = x;
            }
        }
    }
    return vec2();
}
