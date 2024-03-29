#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

vec2 gridpp::downscale_probability(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec2& threshold, const ComparisonOperator& comparison_operator) {
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();
    int nEns = ivalues[0][0].size();

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
            int total = 0;
            for(int k = 0; k < nEns; k++){
                if (gridpp::is_valid(ivalues[I][J][k])) {
                    count = count + 1;
                    if (comparison_operator == gridpp::Leq) {
                        if (ivalues[I][J][k] <= threshold[i][j]) {
                            total = total + 1;
                        }
                    } else if (comparison_operator == gridpp::Lt) {
                        if (ivalues[I][J][k] < threshold[i][j]) {
                            total = total + 1;
                        }
                    } else if (comparison_operator == gridpp::Geq) {
                        if (ivalues[I][J][k] >= threshold[i][j]) {
                            total = total + 1;
                        }
                    } else if (comparison_operator == gridpp::Gt) {
                        if (ivalues[I][J][k] > threshold[i][j]) {
                            total = total + 1;
                        }
                    }   
                }       
            }
            // std::cout << "total[" << i <<"][" << j << "]:" << total << std::endl;
            // std::cout << "count[" << i <<"][" << j << "]:" << count << std::endl;

            float prob;
            if (count == 0) {
                prob = gridpp::MV;
            }
            else {
                prob = (float)total / (float)count;
                // std::cout << "prob[" << i <<"][" << j << "]:" << prob << std::endl;
            }

            output[i][j] = prob;
        }
    }
    return output;
}