#include "gridpp.h"

using namespace gridpp;

vec gridpp::quantile_mapping_curve(const vec& ref, const vec& fcst, vec& output_fcst, vec quantiles) {
    if(ref.size() != fcst.size())
        throw std::invalid_argument("ref and fcst must be of the same size");

    if(quantiles.size() > 0) {
        for(int i = 0; i < quantiles.size(); i++) {
            float curr = quantiles[i];
            if(!gridpp::is_valid(curr) || curr > 1 || curr < 0)
                throw std::invalid_argument("Quantiles must be >= 0 and <= 1");
        }
    }
    output_fcst.clear();
    if(ref.size() == 0) {
        output_fcst = fcst;
        return ref;
    }
    else if(ref.size() == 1) {
        output_fcst = fcst;
        return ref;
    }

    vec ref_sort = ref;
    std::sort(ref_sort.begin(), ref_sort.end());
    vec fcst_sort = fcst;
    std::sort(fcst_sort.begin(), fcst_sort.end());
    int N = quantiles.size();
    int S = fcst_sort.size();
    if(N == 0) {
        output_fcst = fcst_sort;
        return ref_sort;
    }
    else {
        output_fcst.resize(N);
        vec output_ref(N);
        for(int i = 0; i < N; i++) {
            int index = quantiles[i] * (S - 1);
            output_fcst[i] = fcst[index];
            output_ref[i] = ref[index];
        }
        return output_ref;
    }
}
