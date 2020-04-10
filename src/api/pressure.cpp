#include "gridpp.h"
vec2 gridpp::pressure(const gridpp::Grid& igrid, const gridpp::Grid& ogrid, const vec2& ipressure, const vec2& itemperature) {
    vec2 ialtitudes = igrid.get_elevs();
    vec2 oaltitudes = ogrid.get_elevs();
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    int Y = oaltitudes.size();
    vec2 values(Y);
    if(Y > 0) {
        int X = oaltitudes[0].size();
        for(int y = 0; y < Y; y++) {
            values[0].resize(Y, gridpp::MV);
            for(int x = 0; x < X; x++) {
                float temperature = 288.15;
                ivec indices = igrid.get_nearest_neighbour(olats[y][x], olons[y][x]);
                float yy = indices[0];
                float xx = indices[0];
                if(itemperature.size() > 0) {
                    temperature = itemperature[yy][xx];
                }
                values[y][x] = gridpp::pressure(ialtitudes[yy][xx], oaltitudes[y][x], ipressure[y][x], temperature);
            }
        }
    }
    return values;
}
float gridpp::pressure(float ialtitude, float oaltitude, float ipressure, float itemperature) {
    float g0 = 9.80665;
    float M = 0.0289644;
    float R = 8.3144598;
    float value = gridpp::MV;
    if(gridpp::util::is_valid(ialtitude) && gridpp::util::is_valid(oaltitude) && gridpp::util::is_valid(ipressure) && gridpp::util::is_valid(itemperature))
        value = ipressure * exp(-g0 * M * (oaltitude-ialtitude) / (R * itemperature));
    return value;
}

