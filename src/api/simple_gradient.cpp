#include "gridpp.h"

using namespace gridpp;

vec gridpp::simple_gradient(const Grid& igrid, const Points& opoints, const vec2& ivalues, float elev_gradient, Downscaler downscaler) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");

    vec oelevs = opoints.get_elevs();
    vec dvalues = gridpp::downscaling(igrid, opoints, ivalues, downscaler);
    vec delevs = gridpp::downscaling(igrid, opoints, igrid.get_elevs(), downscaler);

    vec output(opoints.size(), gridpp::MV);
    // #pragma omp parallel for
    for(int i = 0; i < opoints.size(); i++) {
        float elev_diff = oelevs[i] - delevs[i];
        float elev_corr = elev_diff * elev_gradient;
        output[i] = dvalues[i] + elev_corr;
    }
    return output;
}

vec2 gridpp::simple_gradient(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, float elev_gradient, Downscaler downscaler) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");

    int nLat = ogrid.size()[0];
    int nLon = ogrid.size()[1];

    vec2 oelevs = ogrid.get_elevs();
    vec2 output = gridpp::init_vec2(nLat, nLon, gridpp::MV);
    vec2 dvalues = gridpp::downscaling(igrid, ogrid, ivalues, downscaler);
    vec2 delevs = gridpp::downscaling(igrid, ogrid, igrid.get_elevs(), downscaler);

    // #pragma omp parallel for
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            float elev_diff = oelevs[i][j] - delevs[i][j];
            float elev_corr = elev_diff * elev_gradient;
            output[i][j] = dvalues[i][j] + elev_corr;
        }
    }
    return output;
}

vec2 gridpp::simple_gradient(const Grid& igrid, const Points& opoints, const vec3& ivalues, float elev_gradient, Downscaler downscaler) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");

    int nPoints = opoints.size();
    int nTime = ivalues.size();

    vec oelevs = opoints.get_elevs();
    vec2 output = gridpp::init_vec2(nTime, nPoints, gridpp::MV);
    vec2 dvalues = gridpp::downscaling(igrid, opoints, ivalues, downscaler);
    vec delevs = gridpp::downscaling(igrid, opoints, igrid.get_elevs(), downscaler);

    #pragma omp parallel for
    for(int i = 0; i < opoints.size(); i++) {
        float elev_diff = oelevs[i] - delevs[i];
        float elev_corr = elev_diff * elev_gradient;
        for(int t = 0; t < nTime; t++) {
            output[t][i] = dvalues[t][i] + elev_corr;
        }
    }
    return output;
}

vec3 gridpp::simple_gradient(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, float elev_gradient, Downscaler downscaler) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");

    int nTime = ivalues.size();
    int nLat = ogrid.size()[0];
    int nLon = ogrid.size()[1];

    vec2 oelevs = ogrid.get_elevs();
    vec3 output = gridpp::init_vec3(nTime, nLat, nLon, gridpp::MV);
    vec3 dvalues = gridpp::downscaling(igrid, ogrid, ivalues, downscaler);
    vec2 delevs = gridpp::downscaling(igrid, ogrid, igrid.get_elevs(), downscaler);

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            float elev_diff = oelevs[i][j] - delevs[i][j];
            float elev_corr = elev_diff * elev_gradient;
            for(int t = 0; t < nTime; t++) {
                output[t][i][j] = dvalues[t][i][j] + elev_corr;
            }
        }
    }
    return output;
}
