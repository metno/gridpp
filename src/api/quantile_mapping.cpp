#include "gridpp.h"

using namespace gridpp;

vec2 gridpp::quantile_mapping_curve(const vec& ref, const vec& fcst, vec quantiles) {
    vec2 curve(2);
    vec ref_sort = ref;
    std::sort(ref_sort.begin(), ref_sort.end());
    vec fcst_sort = fcst;
    std::sort(fcst_sort.begin(), fcst_sort.end());
    int N = quantiles.size();
    int S = fcst_sort.size();
    if(N == 0) {
        curve[0] = fcst_sort;
        curve[1] = ref_sort;
    }
    else {
        curve[0].resize(N);
        curve[1].resize(N);
        for(int i = 0; i < N; i++) {
            int index = quantiles[i] * (S - 1);
            curve[0][i] = fcst[index];
            curve[1][i] = ref[index];
        }
    }
    return curve;
}
