#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

vec2 gridpp::downscale_probability(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec2& threshold, const ComparisonOperator& comparisonoperator) {
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();
    int nEns = ivalues.size();

    vec2 output(nLat);
    for(int i = 0; i < nLat; i++)
        output[i].resize(nLon);

    // std::cout << ivalues.size() << std::endl;

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec indices = igrid.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            int I = indices[0];
            int J = indices[1];
            int count = 0;
            for(int k = 0; k < nEns; k++){
                if (comparisonoperator == gridpp::leq) {
                    if (ivalues[k][I][J] <= threshold[i][j]) {
                        count = count + 1;
                    }
                } else if (comparisonoperator == gridpp::lt) {
                    if (ivalues[k][I][J] < threshold[i][j]) {
                        count = count + 1;
                    }
                } else if (comparisonoperator == gridpp::geq) {
                    if (ivalues[k][I][J] >= threshold[i][j]) {
                        count = count + 1;
                    }
                } else if (comparisonoperator == gridpp::gt) {
                    if (ivalues[k][I][J] > threshold[i][j]) {
                        count = count + 1;
                    }
                }
            // std::cout << "ivalues[" << k << "][" << I << "][" << J << "]:" << ivalues[k][I][J] << std::endl;
            // std::cout << "threshold[" << i << "][" << j << "]:" << threshold[i][j] << std::endl;                
            }
            float prob = (float)count / (float)nEns;
            output[i][j] = prob;
        }
    }
    return output;
}