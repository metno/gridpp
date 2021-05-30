#include "gridpp.h"
#include <math.h>
#include <algorithm>
#include <armadillo>
#include <assert.h>
#include <exception>
#include <boost/math/distributions/normal.hpp>

using namespace gridpp;

namespace {
    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };
}
vec2 gridpp::local_distribution_correction(const Grid& bgrid,
        const vec2& background,
        const Points& points,
        const vec& pobs,
        const vec& pbackground,
        const StructureFunction& structure,
        int max_points) {

    int nY = background.size();
    int nX = background[0].size();
    vec2 blats = bgrid.get_lats();
    vec2 blons = bgrid.get_lons();
    vec2 output(nY);

    float min_quantile = 0.2;
    float max_quantile = 1;
    float min_num_stations = 10;
    bool weighted = true;
    float alpha = 0.01;

    float localizationRadius = structure.localization_distance();

    for(int y = 0; y < nY; y++) {
        output[y].resize(nX, gridpp::MV);
        for(int x = 0; x < nX; x++) {
            output[y][x] = background[y][x];
            if(!gridpp::is_valid(background[y][x]))
                continue;
            float lat = blats[y][x];
            float lon = blons[y][x];
            Point p1 = bgrid.get_point(y, x);
            ivec indices = points.get_neighbours(lat, lon, localizationRadius);

            // Compute fraction of background below current background
            float sum_rho = 0;
            float sum_below = 0;
            int S = indices.size();
            std::vector<std::pair<float, float> > ref_rho;
            std::vector<std::pair<float, float> > fcst_rho;
            ref_rho.reserve(S);
            fcst_rho.reserve(S);
            double sum_ref_rho = 0;
            double sum_fcst_rho = 0;
            int count = 0;
            for(int s = 0; s < S; s++) {
                int index = indices[s];
                if(!gridpp::is_valid(pobs[index]) || !gridpp::is_valid(pbackground[index]))
                    continue;
                if(pobs[index] < 0 || pbackground[index] < 0)
                    continue;
                Point p2 = points.get_point(index);
                float curr_rho = structure.corr_background(p1, p2);
                // curr_rho = 1;
                sum_rho += curr_rho;
                if (pbackground[index] < background[y][x])
                    sum_below += curr_rho;
                float curr_ref = pobs[index];
                float curr_fcst = pbackground[index];
                assert(curr_ref >= 0);
                assert(curr_fcst >= 0);
                ref_rho.push_back(std::pair<float, float>(curr_ref, curr_rho));
                fcst_rho.push_back(std::pair<float, float>(curr_fcst, curr_rho));
                sum_ref_rho += pobs[index] * curr_rho;
                sum_fcst_rho += pbackground[index] * curr_rho;
                count++;
            }
            if (count >= min_num_stations) {
                std::sort(ref_rho.begin(), ref_rho.end(), ::sort_pair_first<float, float>());
                std::sort(fcst_rho.begin(), fcst_rho.end(), ::sort_pair_first<float, float>());

                // Remove the lowest values
                int d0 = (int) count * min_quantile;
                int d1 = count - (int) count * max_quantile;
                ref_rho = std::vector<std::pair<float, float> >(ref_rho.begin() + d0, ref_rho.end() - d1);
                fcst_rho = std::vector<std::pair<float, float> >(fcst_rho.begin() + d0, fcst_rho.end() - d1);

                int new_count = count - d0 + 1 - d1;
                vec ref(new_count, 0);
                vec ref_quantiles(new_count, 0);
                vec fcst(new_count, 0);
                vec fcst_quantiles(new_count, 0);

                assert(ref_rho.size() == new_count - 1);

                // Overwrite the first one to 0
                for(int s = 0; s < new_count - 1; s++) {
                    assert(s < ref_rho.size());
                    ref[s + 1] = ref_rho[s].first;
                    ref_quantiles[s + 1] = ref_quantiles[s] + ref_rho[s].second;
                    fcst[s + 1] = fcst_rho[s].first;
                    fcst_quantiles[s + 1] = fcst_quantiles[s] + fcst_rho[s].second;
                }
                float sum_ref_quantile = ref_quantiles[new_count - 1];
                float sum_fcst_quantile = fcst_quantiles[new_count - 1];

                // Normalize quantiles
                for(int s = 1; s < new_count; s++) {
                    ref_quantiles[s] = min_quantile + ref_quantiles[s] / (sum_ref_quantile) * (max_quantile - min_quantile);
                    fcst_quantiles[s] = min_quantile + fcst_quantiles[s] / (sum_fcst_quantile) * (max_quantile - min_quantile);
                }
                if(background[y][x] > fcst[new_count - 1]) {
                    float diff = ref[new_count - 1] - fcst[new_count - 1];
                    float new_ref = background[y][x] + diff;
                    output[y][x] = new_ref;
                    // output[y][x] = 1;
                }
                else {
                    float q = gridpp::interpolate(background[y][x], fcst, fcst_quantiles);
                    float new_ref = gridpp::interpolate(q, ref_quantiles, ref);
                    // output[y][x] = q;
                    if(weighted) {
                        float w0 = 1 - exp(-alpha * sum_rho);
                        float w1 = 1 - w0;
                        output[y][x] = w0 * new_ref + w1 * background[y][x];
                    }
                    else {
                        output[y][x] = new_ref;
                    }
                }
                // output[y][x] = q;
                // output[y][x] = frac_below * sum_ref_rho;
                // output[y][x] = frac_below; // ref_rho[s].first;
                // output[y][x] = sum_ref_rho / count; // ref_rho[s].first;
                // output[y][x] = acc / count; // ref_rho[s].first;
            }
        }
    }
    return output;
}
