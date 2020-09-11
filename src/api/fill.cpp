#include "gridpp.h"
#include <iostream>

vec2 gridpp::fill(const Grid& igrid, const vec2& input, const Points& points, const vec& radii, float value, bool outside) {
    double s_time = gridpp::util::clock();
    vec lats = points.get_lats();
    vec lons = points.get_lons();
    vec2 output;
    if(outside) {
        output.resize(input.size());
        for(int i = 0; i < input.size(); i++) {
            output[i].resize(input[i].size(), value);
        }
    }
    else {
        output = input;
    }
    for(int i = 0; i < points.size(); i++) {
        ivec2 I = igrid.get_neighbours(lats[i], lons[i], radii[i]);
        // std::cout << " " << i << " " << I.size() << " " << lats[i] << " " << lons[i] << " " << radii[i] << std::endl;
        for(int j = 0; j < I.size(); j++) {
            if(outside)
                output[I[j][0]][I[j][1]] = input[I[j][0]][I[j][1]];
            else
                output[I[j][0]][I[j][1]] = value;
        }
    }
    // std::cout << "Fill time: " << gridpp::util::clock() - s_time << std::endl;
    return output;
}

vec2 gridpp::fill_missing(const vec2& values) {
    int Y = values.size();
    int X = values[0].size();
    vec2 results(Y);
    int count = 0;
    for(int y = 0; y < Y; y++) {
        results[y].resize(X, gridpp::MV);
        int last = 0;
        int next = -1;
        for(int x = 0; x < X; x++) {
            float curr = values[y][x];
            if(!gridpp::util::is_valid(curr)) {
                // Find next
                if(next == -1) {
                    for(next = x; next < X; next++) {
                        if(gridpp::util::is_valid(values[y][next])) {
                            break;
                        }
                    }
                }
                else if(next < x)
                    next = -1;
                // Interpolate
                float value_last = values[y][last];
                float value_next = values[y][next];
                results[y][x] = (value_last) + (value_next - value_last) * (x - last) / (next - last);
                count++;
            }
            else {
                last = x;
                results[y][x] = values[y][x];
            }
        }
    }
    std::cout << "Number of infilled values: " << count << std::endl;
    return results;
}
