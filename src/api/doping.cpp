#include "gridpp.h"

using namespace gridpp;

vec2 gridpp::doping_square(const Grid& igrid, const vec2& background, const Points& points, const vec& observations, const ivec& half_widths, float max_elev_diff) {
    if(gridpp::compatible_size(igrid, background))
        throw std::invalid_argument("Grid size is not the same as observations");
    if(points.size() != observations.size())
        throw std::invalid_argument("Points size is not the same as observations size");
    if(points.size() != half_widths.size())
        throw std::invalid_argument("Points size is not the same as half_widths size");
    if(gridpp::is_valid(max_elev_diff) && max_elev_diff < 0)
        throw std::invalid_argument("max_elev_diff must be greater than or equal to 0");

    double s_time = gridpp::clock();
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();
    const vec2& ielevs = igrid.get_elevs();

    int N = points.size();
    int Y = igrid.size()[0];
    int X = igrid.size()[1];

    for(int i = 0; i < N; i++) {
        if(half_widths[i] < 0)
            throw std::invalid_argument("All half_widths must be greater than or equal to 0");
    }

    vec2 output = background;
    bool check_elev = gridpp::is_valid(max_elev_diff);
    for(int i = 0; i < N; i++) {
        ivec index = igrid.get_nearest_neighbour(lats[i], lons[i]);
        float curr = observations[i];
        for(int yy = std::max(0, index[0] - half_widths[i]); yy <= std::min(Y - 1, index[0] + half_widths[i]); yy++) {
            for(int xx = std::max(0, index[1] - half_widths[i]); xx <= std::min(X - 1, index[0] + half_widths[i]); xx++) {
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

vec2 gridpp::doping_circle(const Grid& igrid, const vec2& background, const Points& points, const vec& observations, const vec& radii, float max_elev_diff) {
    if(gridpp::compatible_size(igrid, background))
        throw std::invalid_argument("Grid size is not the same as observations");
    if(points.size() != observations.size())
        throw std::invalid_argument("Points size is not the same as observations size");
    if(points.size() != radii.size())
        throw std::invalid_argument("Points size is not the same as radii size");
    if(gridpp::is_valid(max_elev_diff) && max_elev_diff < 0)
        throw std::invalid_argument("max_elev_diff must be greater than or equal to 0");

    double s_time = gridpp::clock();
    const vec& lats = points.get_lats();
    const vec& lons = points.get_lons();
    const vec& elevs = points.get_elevs();
    const vec2& ielevs = igrid.get_elevs();

    int N = points.size();
    int Y = igrid.size()[0];
    int X = igrid.size()[1];

    for(int i = 0; i < N; i++) {
        if(radii[i] < 0)
            throw std::invalid_argument("radii must be greater than or equal to 0");
    }

    vec2 output = background;
    bool check_elev = gridpp::is_valid(max_elev_diff);
    for(int i = 0; i < N; i++) {
        ivec2 I = igrid.get_neighbours(lats[i], lons[i], radii[i]);
        float curr = observations[i];
        for(int j = 0; j < I.size(); j++) {
            int yy = I[j][0];
            int xx = I[j][1];
            if(check_elev) {
                float diff = fabs(elevs[i] - ielevs[yy][xx]);
                if(diff > max_elev_diff)
                    continue;
            }
            output[yy][xx] = curr;
        }
    }
    return output;

}
