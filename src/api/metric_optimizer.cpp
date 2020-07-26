#include "gridpp.h"
#include <iostream>
#include <algorithm>
#include <boost/math/tools/minima.hpp>
using boost::math::tools::brent_find_minima;

namespace {
    float partial(float ref_diff, float fcst_diff, float sigma) {
        float total = 0;
        float ref_total = 0;
        float fcst_total = 1;
        if(ref_diff > 0)
            ref_total = 1;
        else
            ref_total = 1 + (ref_diff) / sigma;
        if(fcst_diff > 0)
            fcst_total = 1;
        else
            fcst_total = 1 + (fcst_diff) / sigma;

        return ref_total * fcst_total;
    }
    class Data {
        public:
            Data(const vec& ref, const vec& fcst, float threshold, gridpp::Metric metric) {
                this->ref = ref;
                this->fcst = fcst;
                this->metric = metric;
                this->threshold = threshold;
                int N = ref.size();
                min = gridpp::MV;
                max = gridpp::MV;
                for(int i = 0; i < N; i++) {
                    if(!gridpp::util::is_valid(min) || fcst[i] < min) {
                        min = fcst[i];
                    }
                    if(!gridpp::util::is_valid(max) || fcst[i] > max) {
                        max = fcst[i];
                    }
                }
            };
            vec ref;
            vec fcst;
            float threshold;
            gridpp::Metric metric;
            float min;
            float max;
    };
    class Score {
        public:
            Score(float a, float b, float c, float d, gridpp::Metric metric) {
                this->a = a;
                this->b = b;
                this->c = c;
                this->d = d;
                this->metric = metric;
            };
            float a;
            float b;
            float c;
            float d;
            gridpp::Metric metric;
    };
    class Score_func {
        public:
            Score_func(const Data& idata) : data(idata) {
            };
            float operator()(double threshold) {
                float a = 0, b = 0, c = 0, d = 0;
                int N = data.fcst.size();
                for(int t = 0; t < N; t++) {
                    if((data.ref[t] > data.threshold) && (data.fcst[t] > threshold))
                        a++;
                    else if((data.ref[t] <= data.threshold) && (data.fcst[t] > threshold))
                        b++;
                    else if((data.ref[t] > data.threshold) && (data.fcst[t] <= threshold))
                        c++;
                    else if((data.ref[t] <= data.threshold) && (data.fcst[t] <= threshold))
                        d++;
                }
                float s = - gridpp::calc_score(a, b, c, d, data.metric);
                if(s > -0.0001)
                    s = 0;

                // std::cout << "Testing: " << threshold << " = " << s << std::endl;
                return s;
            }
            float calc_uncertainty(double threshold) {
                float a = 0, b = 0, c = 0, d = 0;
                int N = data.fcst.size();
                for(int t = 0; t < N; t++) {
                    if((data.ref[t] > data.threshold) && (data.fcst[t] > threshold))
                        a++;
                    else if((data.ref[t] <= data.threshold) && (data.fcst[t] > threshold))
                        b++;
                    else if((data.ref[t] > data.threshold) && (data.fcst[t] <= threshold))
                        c++;
                    else if((data.ref[t] <= data.threshold) && (data.fcst[t] <= threshold))
                        d++;
                }
                float s = - gridpp::calc_score(a, b, c, d, data.metric);

                float zero = 0;
                float s0 = - gridpp::calc_score(std::max(zero, a-1), std::max(zero, b-1), std::max(zero, c-1), std::max(zero, d-1), data.metric);
                float diff = fabs(fabs(s - s0) / std::min(s, s0));
                if(diff == -gridpp::MV)
                    diff = gridpp::MV;
                return diff;
            }
            const Data& data;
    };
}

vec2 gridpp::metric_optimizer_curve(const vec& ref, const vec& fcst, const vec& thresholds, gridpp::Metric metric) {
    int N = thresholds.size();
    vec x(N);
    vec y(N);
    float min = 0;
    float max = 10;
    vec2 results(2);
    results[1] = thresholds;
    results[0].resize(N);

    float last_diff = 0;
    for(int i = 0; i < N; i++) {
        // Compute score for different combinations
        // vec[i] = thresholds[i];
        vec scores;
        float value = gridpp::get_optimal_threshold(ref, fcst, thresholds[i], metric, scores);
        results[0][i] = value;
        last_diff = value - thresholds[i];
    }
    return results;
}

vec gridpp::compute_scores(const vec& ref, const vec& fcst, float othreshold, const vec& fthresholds, float sigma, gridpp::Metric metric) {
    int num_tests = fthresholds.size();
    vec scores(num_tests, gridpp::MV);
    int N = ref.size();
    float threshold = othreshold;
    for(int j = 0; j < num_tests; j++) {
        float fthreshold = fthresholds[j];
        float a = 0, b = 0, c = 0, d = 0;
        for(int t = 0; t < N; t++) {
            if(sigma == 0) {
                if((ref[t] > threshold) && (fcst[t] > fthreshold))
                    a++;
                else if((ref[t] <= threshold) && (fcst[t] > fthreshold))
                    b++;
                else if((ref[t] > threshold) && (fcst[t] <= fthreshold))
                    c++;
                else if((ref[t] <= threshold) && (fcst[t] <= fthreshold))
                    d++;
            }
            else {
                float ref_diff = ref[t] - threshold;
                float fcst_diff = fcst[t] - fthreshold;
                if((ref[t] > threshold - sigma) && (fcst[t] > fthreshold - sigma))
                    a += ::partial(ref_diff, fcst_diff, sigma);
                if((ref[t] <= threshold + sigma) && (fcst[t] > fthreshold - sigma))
                    b += ::partial(-ref_diff, fcst_diff, sigma);
                if((ref[t] > threshold - sigma) && (fcst[t] <= fthreshold + sigma))
                    c += ::partial(ref_diff, -fcst_diff, sigma);
                if((ref[t] <= threshold + sigma) && (fcst[t] <= fthreshold + sigma))
                    d += ::partial(-ref_diff, -fcst_diff, sigma);

            }
        }
        if(0 && (a < 2 || b < 2 || c < 2 || d < 2))
            scores[j] = gridpp::MV;
        else {
            float s = calc_score(a, b, c, d, metric);
            float zero = 0;
            float s0 = calc_score(std::max(zero, a-1), std::max(zero, b-1), std::max(zero, c-1), std::max(zero, d-1), metric);
            // float s1 = calc_score(a-1, b-1, c-1, d-1, metric);
            float d = fabs(s - s0) / std::max(s, s0);
            // std::cout << s << " " << s0 << " " << s1 << " " << fabs(s - s0) / s << std::endl;
            if(!gridpp::util::is_valid(s) || !gridpp::util::is_valid(s0))
                scores[j] = gridpp::MV;
            else if(gridpp::util::is_valid(d) && d > 1e-3)
                scores[j] = gridpp::MV;
            // if(a < 20 || b < 20 || c < 20 || d < 20)
            //     scores[j] = gridpp::MV;
            // else
            else if(s <= 0.0001)
                scores[j] = 0;
            else
                scores[j] = s;
        }
    }
    return scores;
}

float gridpp::get_optimal_threshold2(const vec& ref, const vec& fcst, float threshold, gridpp::Metric metric) {
    // int bits = std::numeric_limits<double>::digits;
    int bits = 32;
    Data data(ref, fcst, threshold, metric);
    Score_func func(data);
    std::pair<double, double> r = brent_find_minima<Score_func, double>(func, data.min, data.max, bits);
    float x = r.first;
    if(r.second >= -0.0001)
        x = gridpp::MV;
    float diff = func.calc_uncertainty(r.first);
    float s0 = func(data.min);
    float s1 = func(data.max);
    // Remove points where uncertainty is large
    return x;
    if(!gridpp::util::is_valid(r.second))
        x = gridpp::MV;
    else if(r.second >= -0.0001)
        x = gridpp::MV;
    else if(!gridpp::util::is_valid(diff))
        x = gridpp::MV;
    else if(diff > 1e-3)
        x = gridpp::MV;
    else if (fabs(r.second - s0) < 0.001 || fabs(r.second - s1) < 0.001)
        x = gridpp::MV;
    std::cout << r.first << " " << r.second << " " << diff << " " << x << " " << s0 << " " << s1 << std::endl;
    return x;
}
float gridpp::get_optimal_threshold(const vec& ref, const vec& fcst, float threshold, gridpp::Metric metric, vec& scores) {
    int N = ref.size();
    float min_ref = gridpp::MV;
    float max_ref = gridpp::MV;
    for(int i = 0; i < N; i++) {
        if(!gridpp::util::is_valid(min_ref) || fcst[i] < min_ref) {
            min_ref = fcst[i];
        }
        if(!gridpp::util::is_valid(max_ref) || fcst[i] > max_ref) {
            max_ref = fcst[i];
        }
    }

    int num_tests = 10;
    vec bins(num_tests);
    float step = 0.1;
    for(int j = 0; j < num_tests; j++) {
        bins[j] = (j - num_tests / 2.0) * step;
        // std::cout << bins[j] << std::endl;
    }
    scores.clear();
    scores.resize(num_tests);
    float optimal = gridpp::MV;
    float max_score = gridpp::MV;
    float min_score = gridpp::MV;
    for(int j = 0; j < num_tests; j++) {
        float fthreshold = threshold + bins[j];
        if(fthreshold < min_ref || fthreshold > max_ref)
            continue;
        float a = 0, b = 0, c = 0, d = 0;
        for(int t = 0; t < N; t++) {
            if(0) {
                if((ref[t] > threshold) && (fcst[t] > fthreshold))
                    a++;
                else if((ref[t] <= threshold) && (fcst[t] > fthreshold))
                    b++;
                else if((ref[t] > threshold) && (fcst[t] <= fthreshold))
                    c++;
                else if((ref[t] <= threshold) && (fcst[t] <= fthreshold))
                    d++;
            }
            else {
                float sigma = 0.05;
                float ref_diff = ref[t] - threshold;
                float fcst_diff = fcst[t] - fthreshold;
                if((ref[t] > threshold - sigma) && (fcst[t] > fthreshold - sigma))
                    a += ::partial(ref_diff, fcst_diff, sigma);
                if((ref[t] <= threshold + sigma) && (fcst[t] > fthreshold - sigma))
                    b += ::partial(-ref_diff, fcst_diff, sigma);
                if((ref[t] > threshold - sigma) && (fcst[t] <= fthreshold + sigma))
                    c += ::partial(ref_diff, -fcst_diff, sigma);
                if((ref[t] <= threshold + sigma) && (fcst[t] <= fthreshold + sigma))
                    d += ::partial(-ref_diff, -fcst_diff, sigma);

            }
        }
        std::cout << a << " " << b << " " << c << " " << d << std::endl;
        // if(a + b < 50 || c + d < 50)
        //     continue;

        float s = calc_score(a, b, c, d, metric);
        scores[j] = s;
        std::cout << fthreshold << " " << s << std::endl;
        if(!gridpp::util::is_valid(max_score) || (s >= max_score)) {
            max_score = s;
            optimal = fthreshold;
        }
        if(!gridpp::util::is_valid(min_score) || (s < min_score)) {
            min_score = s;
        }
        // std::cout << fthreshold << " " << s << " " << max_score << " " << min_score << " " << a << " " << b << " " << c << " " << d << std::endl;
    }
    if(0) {

    }
    else if(optimal == threshold + bins[0] || optimal == threshold + bins[bins.size() - 1])
        optimal = gridpp::MV;
    else if(!gridpp::util::is_valid(max_score))
        optimal = gridpp::MV;
    else if(gridpp::util::is_valid(min_score) && gridpp::util::is_valid(max_score) && (max_score - min_score < 0.01))
        optimal = gridpp::MV;
    else if(max_score - scores[0] < 0.001)
        optimal = gridpp::MV;
    else if(max_score < 0.0001)
        optimal = gridpp::MV;
    // std::cout << threshold << " " << optimal << " " << max_score << std::endl;
    return optimal;
}
float gridpp::calc_score(float a, float b, float c, float d, gridpp::Metric metric) {
    if(metric == gridpp::Ets) {
        float N = a + b + c + d;
        float ar = (a + b) / 1.0 / N * (a + c);
        if(a + b + c - ar == 0)
            return gridpp::MV;
        return (a - ar) / 1.0 / (a + b + c - ar);
    }
    else if(metric == gridpp::Ts) {
        float N = a + b + c + d;
        return a / 1.0 / (a + b + c);
    }
    else if(metric == gridpp::Pc) {
        float N = a + b + c + d;
        return (a + d) / N;
    }
    else if(metric == gridpp::Kss) {
        if((a + c) * (b + d) == 0)
            return gridpp::MV;
        return (a * d - b * c) * 1.0 / ((a + c) * (b + d));
    }
    else if(metric == gridpp::Bias) {
        if(b + c == 0)
            return gridpp::MV;
        float value = -float(fabs(b - c))/(b + c);
        return value;
    }
    else if(metric == gridpp::Hss) {
        float denom = ((a + c) * (c + d) + (a + b) * (b + d));
        if(denom == 0)
            return gridpp::MV;
        return 2.0 * (a * d - b * c) / denom;
    }
}
ivec gridpp::monotonize(const vec2& curve) {
    vec2 new_curve;
    ivec new_indices;
    assert(curve.size() == 2);
    int N = curve[0].size();
    int Nmiddle = N / 2;
    new_curve.reserve(N);
    new_indices.reserve(N);
    // Upper half
    int start = 0;
    float last = curve[0][start];
    bool deviation = false;
    int deviation_index = 0;
    float x_min = 0;
    float x_max = 0;
    float x_min_index = 0;
    for(int i = start; i < N; i++) {
        float x = curve[0][i];
        float y = curve[1][i];
        if(!gridpp::util::is_valid(x) || !gridpp::util::is_valid(y))
            continue;
        if(deviation) {
            if(x < x_min) {
                x_min = x;
                x_min_index = i;
            }
            if(x > x_max) {
                // Stopping criterion
                x_max = x;
                // std::cout << "End of loop at " << i << " x=" << x << "(" << x_min << "," << x_max << ")" << std::endl;
                // Find latest x point
                for(int j = new_indices.size() - 1; j >= 0; j--) {
                    int index = new_indices[j];
                    if(curve[0][index] < x_min)
                        break;
                    else {
                        // std::cout << "Removing index=" << index << " x=" << curve[0][index] << std::endl;
                        new_indices.pop_back();
                    }
                }
                new_indices.push_back(i);
            }
            deviation = false;
        }
        else {
            if(x <= last) {
                // We have a deviation
                deviation = true;
                deviation_index = i;
                x_min = x;
                x_max = x;
                x_min_index = i;
                // std::cout << "Detecting loop from " << i << " x=" << x << std::endl;
            }
            else {
                new_indices.push_back(i);
                last = x;
            }
        }
    }
    return new_indices;
}
