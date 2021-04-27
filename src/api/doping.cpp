#include "gridpp.h"

using namespace gridpp;

vec2 gridpp::doping(const Grid& igrid, const vec2& input, const Points& points, const vec& values, int half_width, float max_elev_diff) {
    if(gridpp::compatible_size(igrid, input))
        throw std::invalid_argument("Grid size is not the same as values");
    if(points.size() != values.size())
        throw std::invalid_argument("Points size is not the same as values size");
    if(gridpp::is_valid(max_elev_diff) && max_elev_diff < 0)
        throw std::invalid_argument("max_elev_diff must be greater than or equal to 0");
    if(half_width < 0)
        throw std::invalid_argument("half_width must be greater than or equal to 0");

    double s_time = gridpp::clock();
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();
    const vec2& ielevs = igrid.get_elevs();

    int N = points.size();
    int Y = igrid.size()[0];
    int X = igrid.size()[1];

    vec2 output = input;
    bool check_elev = gridpp::is_valid(max_elev_diff);
    for(int i = 0; i < N; i++) {
        ivec index = igrid.get_nearest_neighbour(lats[i], lons[i]);
        float curr = values[i];
        for(int yy = std::max(0, index[0] - half_width); yy <= std::min(Y - 1, index[0] + half_width); yy++) {
            for(int xx = std::max(0, index[1] - half_width); xx <= std::min(X - 1, index[0] + half_width); xx++) {
                if(check_elev) {
                    float diff = fabs(elevs[i] - ielevs[yy][xx]);
                    if(diff > max_elev_diff)
                        continue;
                }
                output[yy][xx] = curr;
            }
        }
    }
    return output;

}
