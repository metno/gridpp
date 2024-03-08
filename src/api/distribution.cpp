#include "gridpp.h"

using namespace gridpp;

vec gridpp::gamma_inv(const vec& levels, const vec& shape, const vec& scale) {
    int N = shape.size();
    for(int i = 0; i < N; i++) {
        if(levels[i] < 0 || levels[i] > 1 || !gridpp::is_valid(levels[i])) {
            std::stringstream ss;
            ss << "Invalid level '" << levels[i] << "'. Levels must be on the interval [0, 1].";
            throw std::invalid_argument(ss.str());
        }
        if(shape[i] <= 0 || !gridpp::is_valid(shape[i])) {
            std::stringstream ss;
            ss << "Invalid shape '" << shape[i] << "'. Shapes must be > 0.";
            throw std::invalid_argument(ss.str());
        }
        if(scale[i] <= 0 || !gridpp::is_valid(scale[i])) {
            std::stringstream ss;
            ss << "Invalid scale '" << scale[i] << "'. Scale must be > 0.";
            throw std::invalid_argument(ss.str());
        }
    }
    vec results(N);

    #pragma omp parallel for
    for(int i = 0; i < N; i++) {
        boost::math::gamma_distribution<> m_gamma_dist(shape[i], scale[i]);
        float value = boost::math::quantile(m_gamma_dist, levels[i]);
        results[i] = value;
    }
    return results;
}
