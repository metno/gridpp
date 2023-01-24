#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

namespace{
    vec2 mask_threshold_downscale(const Grid& igrid, const Grid& ogrid, const vec3& ivalues_true, const vec3& ivalues_false, const vec3& threshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const Statistic& statistic, const float quantile);
}

vec2 gridpp::mask_threshold_downscale_consensus(const Grid& igrid, const Grid& ogrid, const vec3& ivalues_true, const vec3& ivalues_false, const vec3& threshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const Statistic& statistic) {
    return ::mask_threshold_downscale(igrid, ogrid, ivalues_true, ivalues_false, threshold_values, threshold, comparison_operator, statistic, 0);
}
vec2 gridpp::mask_threshold_downscale_quantile(const Grid& igrid, const Grid& ogrid, const vec3& ivalues_true, const vec3& ivalues_false, const vec3& threshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const float quantile) {
    return ::mask_threshold_downscale(igrid, ogrid, ivalues_true, ivalues_false, threshold_values, threshold, comparison_operator, gridpp::Quantile, quantile);
}

namespace {
    vec2 mask_threshold_downscale(const Grid& igrid, const Grid& ogrid, const vec3& ivalues_true, const vec3& ivalues_false, const vec3& threshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const Statistic& statistic, const float quantile) {
        vec2 iOutputLats = ogrid.get_lats();
        vec2 iOutputLons = ogrid.get_lons();

        int nLat = iOutputLats.size();
        int nLon = iOutputLats[0].size();
        int nEns = ivalues_true[0][0].size();

        vec2 output(nLat);
        for(int i = 0; i < nLat; i++)
            output[i].resize(nLon);

        #pragma omp parallel for collapse(2)
        for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
                ivec indices = igrid.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
                int I = indices[0];
                int J = indices[1];
                int count = 0;
                vec masked_ivalues(nEns, gridpp::MV);
                for(int k = 0; k < nEns; k++) {
                    if (gridpp::is_valid(threshold_values[I][J][k] )) {
                        if (comparison_operator == gridpp::Leq) {
                            if (threshold_values[I][J][k] <= threshold[i][j]) {
                                masked_ivalues[k] = ivalues_true[I][J][k];
                            }
                            else {
                                masked_ivalues[k] = ivalues_false[I][J][k];
                            }
                        } else if (comparison_operator == gridpp::Lt) {
                            if (threshold_values[I][J][k] < threshold[i][j]) {
                                masked_ivalues[k] = ivalues_true[I][J][k];
                            }
                            else {
                                masked_ivalues[k] = ivalues_false[I][J][k];
                            }
                        } else if (comparison_operator == gridpp::Geq) {
                            if (threshold_values[I][J][k] >= threshold[i][j]) {
                                masked_ivalues[k] = ivalues_true[I][J][k];
                            }
                            else {
                                masked_ivalues[k] = ivalues_false[I][J][k];
                            }
                        } else if (comparison_operator == gridpp::Gt) {
                            if (threshold_values[I][J][k] > threshold[i][j]) {
                                masked_ivalues[k] = ivalues_true[I][J][k];
                            }
                            else {
                                masked_ivalues[k] = ivalues_false[I][J][k];
                            }
                        }
                    }

                }
                if(statistic == gridpp::Quantile) {
                    output[i][j] = gridpp::calc_quantile(masked_ivalues, quantile);
                }
                else {
                    output[i][j] = gridpp::calc_statistic(masked_ivalues, statistic);
                }
            }
        }
        return output;
    }
}
