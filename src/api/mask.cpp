#include "gridpp.h"

vec2 gridpp::mask(const Grid& igrid, const vec2& input, const Points& points, const vec& radii, float value, bool keep) {
    double s_time = gridpp::util::clock();
    vec lats = points.get_lats();
    vec lons = points.get_lons();
    vec2 output;
    if(keep) {
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
            if(keep)
                output[I[j][0]][I[j][1]] = input[I[j][0]][I[j][1]];
            else
                output[I[j][0]][I[j][1]] = value;
        }
    }
    std::cout << "Mask time: " << gridpp::util::clock() - s_time << std::endl;
    return output;
}
