#include "gridpp.h"

using namespace gridpp;

vec gridpp::gamma_inv(const vec& levels, const vec& shape, const vec& scale) {
    int N = shape.size();
    vec results(N);

    #pragma omp parallel for
    for(int i = 0; i < N; i++) {
        boost::math::gamma_distribution<> m_gamma_dist(shape[i], scale[i]);
        float value = boost::math::quantile(m_gamma_dist, levels[i]);
        results[i] = value;
    }
    return results;
}
