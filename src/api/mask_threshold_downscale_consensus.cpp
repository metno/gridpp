#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

vec2 gridpp::mask_threshold_downscale_consensus(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec3& threshold_values, const vec2& threshold, const ComparisonOperator& comparisonoperator, const Statistic& statistic) {
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();
    int nEns = ivalues.size();

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
            vec masked_ivalues(nEns);
            for(int k = 0; k < nEns; k++){
                if (gridpp::is_valid(threshold_values[k][I][J] )) {
                    if (comparisonoperator == gridpp::leq) {
                        if (threshold_values[k][I][J] <= threshold[i][j]) {
                            masked_ivalues[k] = ivalues[k][I][J];
                        }
                        else {
                            masked_ivalues[k] = 0.0;
                        }
                    } else if (comparisonoperator == gridpp::lt) {
                        if (threshold_values[k][I][J] < threshold[i][j]) {
                            masked_ivalues[k] = ivalues[k][I][J];
                        }
                        else {
                            masked_ivalues[k] = 0.0;
                        }
                    } else if (comparisonoperator == gridpp::geq) {
                        if (threshold_values[k][I][J] >= threshold[i][j]) {
                            masked_ivalues[k] = ivalues[k][I][J];
                        }
                        else {
                            masked_ivalues[k] = 0.0;
                        }
                    } else if (comparisonoperator == gridpp::gt) {
                        if (threshold_values[k][I][J] > threshold[i][j]) {
                            masked_ivalues[k] = ivalues[k][I][J];
                        }
                        else {
                            masked_ivalues[k] = 0.0;
                        }
                    }    
                } 
                else {
                    masked_ivalues[k] = NAN;
                }
                // std::cout << "masked_ivalues[" << i <<"][" << j << "][" << k << "]" << masked_ivalues[k] << std::endl; 
            }
            output[i][j] = gridpp::calc_statistic(masked_ivalues, statistic);
            // std::cout << "count[" << i <<"][" << j << "]:" << output[i][j] << std::endl;
        }
    }
    return output;
}