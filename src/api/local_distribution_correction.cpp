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
        float min_quantile,
        float max_quantile,
        int min_points) {
    vec2 pobs0;
    pobs0.push_back(pobs);
    vec2 pbackground0;
    pbackground0.push_back(pbackground);
    return gridpp::local_distribution_correction(bgrid, background, points, pobs0, pbackground0, structure, min_quantile, max_quantile, min_points);
}
vec2 gridpp::local_distribution_correction(const Grid& bgrid,
        const vec2& background,
        const Points& points,
        const vec2& pobs,
        const vec2& pbackground,
        const StructureFunction& structure,
        float min_quantile,
        float max_quantile,
        int min_points) {

    int nY = background.size();
    int nX = background[0].size();
    int nT = pobs.size();
    vec2 blats = bgrid.get_lats();
    vec2 blons = bgrid.get_lons();
    vec2 output(nY);

    if (pobs.size() != pbackground.size()) {
        std::stringstream ss;
        ss << "pobs (" << pobs.size() << "," << pobs[0].size() << ") is not the same size as pbackground (" << pbackground.size() << "," << pbackground[0].size() << ")";
        throw std::invalid_argument(ss.str());
    }

    bool weighted = true;
    float alpha = 0.01;

    float localizationRadius = structure.localization_distance();

    for(int y = 0; y < nY; y++) {
        output[y].resize(nX, gridpp::MV);
        for(int x = 0; x < nX; x++) {
            // Default to the background value, if there are no observations
            output[y][x] = background[y][x];

            if(!gridpp::is_valid(background[y][x]))
                continue;

            // FInd all stations within the localization radius of the structure function
            float lat = blats[y][x];
            float lon = blons[y][x];
            Point p1 = bgrid.get_point(y, x);
            ivec indices = points.get_neighbours(lat, lon, localizationRadius);

            float sum_below = 0;
            int nS = indices.size();
            std::vector<std::pair<float, float> > ref_rho;
            std::vector<std::pair<float, float> > fcst_rho;
            ref_rho.reserve(nS * nT);
            fcst_rho.reserve(nS * nT);
            float sum_rho = 0;
            int count = 0;

            // Put ref and background into arrays
            for(int s = 0; s < nS; s++) {
                int index = indices[s];
                Point p2 = points.get_point(index);
                float curr_rho = structure.corr_background(p1, p2);
                // curr_rho = 1;
                for(int t = 0; t < nT; t++) {
                    if(!gridpp::is_valid(pobs[t][index]) || !gridpp::is_valid(pbackground[t][index]))
                        continue;
                    if(pobs[t][index] < 0 || pbackground[t][index] < 0)
                        continue;
                    float curr_ref = pobs[t][index];
                    float curr_fcst = pbackground[t][index];
                    assert(curr_ref >= 0);
                    assert(curr_fcst >= 0);
                    ref_rho.push_back(std::pair<float, float>(curr_ref, curr_rho));
                    fcst_rho.push_back(std::pair<float, float>(curr_fcst, curr_rho));
                    sum_rho += curr_rho;
                    count++;
                }
            }

            if (count >= min_points) {
                // Create a calibration curve of ref,fcst, so that we can adjust the current
                // background value
                std::sort(ref_rho.begin(), ref_rho.end(), ::sort_pair_first<float, float>());
                std::sort(fcst_rho.begin(), fcst_rho.end(), ::sort_pair_first<float, float>());

                // Remove the lowest and highest values
                int d0 = (int) count * min_quantile;
                int d1 = count - (int) count * max_quantile;
                ref_rho = std::vector<std::pair<float, float> >(ref_rho.begin() + d0, ref_rho.end() - d1);
                fcst_rho = std::vector<std::pair<float, float> >(fcst_rho.begin() + d0, fcst_rho.end() - d1);

                // Don't create precip out of thin air
                if(fcst_rho[fcst_rho.size() - 1].first < 0.1) {
                    continue;
                }

                // Create a new calculation curve, add an extra point to 0,0.
                int new_count = count - d0 + 1 - d1;
                vec ref(new_count, 0);
                vec ref_quantiles(new_count, 0);
                vec fcst(new_count, 0);
                vec fcst_quantiles(new_count, 0);

                assert(ref_rho.size() == new_count - 1);

                for(int s = 0; s < new_count - 1; s++) {
                    assert(s < ref_rho.size());
                    ref[s + 1] = ref_rho[s].first;
                    ref_quantiles[s + 1] = ref_quantiles[s] + ref_rho[s].second;
                    fcst[s + 1] = fcst_rho[s].first;
                    fcst_quantiles[s + 1] = fcst_quantiles[s] + fcst_rho[s].second;
                }
                float sum_ref_quantile = ref_quantiles[new_count - 1];
                float sum_fcst_quantile = fcst_quantiles[new_count - 1];

                // Normalize quantiles to be between min_quantile and max_quantile
                for(int s = 1; s < new_count; s++) {
                    ref_quantiles[s] = min_quantile + ref_quantiles[s] / (sum_ref_quantile) * (max_quantile - min_quantile);
                    fcst_quantiles[s] = min_quantile + fcst_quantiles[s] / (sum_fcst_quantile) * (max_quantile - min_quantile);
                }

                // Correct values above the curve by maintaining the bias at the end of the curve
                if(background[y][x] > fcst[new_count - 1]) {
                    float diff = ref[new_count - 1] - fcst[new_count - 1];
                    float new_ref = background[y][x] + diff;
                    output[y][x] = new_ref;
                }
                else {
                    float q = gridpp::interpolate(background[y][x], fcst, fcst_quantiles);
                    float new_ref = gridpp::interpolate(q, ref_quantiles, ref);
                    if(weighted) {
                        float w0 = 1 - exp(-alpha * sum_rho);
                        float w1 = 1 - w0;
                        output[y][x] = w0 * new_ref + w1 * background[y][x];
                    }
                    else {
                        output[y][x] = new_ref;
                    }
                }
            }
        }
    }
    return output;
}
