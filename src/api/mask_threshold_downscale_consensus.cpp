#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

vec2 gridpp::mask_threshold_downscale_consensus(const Grid& igrid, const Grid& ogrid, const vec3& ivaluestrue, const vec3& ivaluesfalse, const vec3& threshold_values, const vec2& threshold, const ComparisonOperator& comparison_operator, const Statistic& statistic) {
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();
    int nEns = ivaluestrue[0][0].size();

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
            for(int k = 0; k < nEns; k++){
                if (gridpp::is_valid(threshold_values[I][J][k] )) {
                    if (comparison_operator == gridpp::Leq) {
                        if (threshold_values[I][J][k] <= threshold[i][j]) {
                            masked_ivalues[k] = ivaluestrue[I][J][k];
                        }
                        else {
                            masked_ivalues[k] = ivaluesfalse[I][J][k];
                        }
                    } else if (comparison_operator == gridpp::Lt) {
                        if (threshold_values[I][J][k] < threshold[i][j]) {
                            masked_ivalues[k] = ivaluestrue[I][J][k];
                        }
                        else {
                            masked_ivalues[k] = ivaluesfalse[I][J][k];
                        }
                    } else if (comparison_operator == gridpp::Geq) {
                        if (threshold_values[I][J][k] >= threshold[i][j]) {
                            masked_ivalues[k] = ivaluestrue[I][J][k];
                        }
                        else {
                            masked_ivalues[k] = ivaluesfalse[I][J][k];
                        }
                    } else if (comparison_operator == gridpp::Gt) {
                        if (threshold_values[I][J][k] > threshold[i][j]) {
                            masked_ivalues[k] = ivaluestrue[I][J][k];
                        }
                        else {
                            masked_ivalues[k] = ivaluesfalse[I][J][k];
                        }
                    }    
                }

            }
            output[i][j] = gridpp::calc_statistic(masked_ivalues, statistic);                
        }
    }
    return output;
}