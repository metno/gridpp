#include "gridpp.h"
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <assert.h>
#include <execinfo.h>
#include <signal.h>
#include <iomanip>
#include <cstdio>
#include <exception>

#ifdef DEBUG
extern "C" void __gcov_flush();
#endif

bool gridpp::util::is_valid(float value) {
    return !std::isnan(value) && !std::isinf(value) && value != gridpp::MV && value != MV_CML;
}
float gridpp::util::calc_statistic(const vec& array, gridpp::Statistic statistic) {
    // Initialize to missing
    float value = gridpp::MV;
    if(statistic == gridpp::Mean || statistic == gridpp::Sum) {
        float total = 0;
        int count = 0;
        for(int n = 0; n < array.size(); n++) {
            if(gridpp::util::is_valid(array[n])) {
                total += array[n];
                count++;
            }
        }
        if(count > 0) {
            if(statistic == gridpp::Mean)
                value = total / count;
            else
                value = total;
        }
    }
    else if(statistic == gridpp::Std || statistic == gridpp::Variance) {
        // STD = sqrt(E[X^2] - E[X]^2)
        // The above formula is unstable when the variance is small and the mean is large.
        // Use the property that VAR(X) = VAR(X-K). Provided K is any element in the array,
        // the resulting calculation of VAR(X-K) is stable. Set K to the first non-missing value.
        float total  = 0;
        float total2 = 0;
        float K = gridpp::MV;
        int count = 0;
        for(int n = 0; n < array.size(); n++) {
            if(gridpp::util::is_valid(array[n])) {
                if(!gridpp::util::is_valid(K))
                    K = array[n];
                assert(gridpp::util::is_valid(K));

                total  += array[n] - K;
                total2 += (array[n] - K)*(array[n] - K);
                count++;
            }
        }
        if(count > 0) {
            float mean  = total / count;
            float mean2 = total2 / count;
            float var   = mean2 - mean*mean;
            if(var < 0) {
                // This should never happen
                var = 0;
                // Util::warning("CalibratorNeighbourhood: Problems computing std, unstable result. Setting value to 0");
            }
            value = var;
            if(statistic == gridpp::Std) {
                value = sqrt(var);
            }
        }
    }
    else {
        float quantile = gridpp::MV;
        if(statistic == gridpp::Min)
            quantile = 0;
        else if(statistic == gridpp::Median)
            quantile = 0.5;
        else if(statistic == gridpp::Max)
            quantile = 1;
        else
            throw std::runtime_error("Internal error. Cannot compute statistic");
        value = gridpp::util::calc_quantile(array, quantile);
    }
    return value;
}
float gridpp::util::calc_quantile(const vec& array, float quantile) {
    // Initialize to missing
    float value = gridpp::MV;
    // Remove missing
    std::vector<float> cleanHood;
    cleanHood.reserve(array.size());
    for(int i = 0; i < array.size(); i++) {
        if(gridpp::util::is_valid(array[i]))
            cleanHood.push_back(array[i]);
    }
    int N = cleanHood.size();
    if(N > 0) {
        std::sort(cleanHood.begin(), cleanHood.end());
        int lowerIndex = floor(quantile * (N-1));
        int upperIndex = ceil(quantile * (N-1));
        float lowerQuantile = (float) lowerIndex / (N-1);
        float upperQuantile = (float) upperIndex / (N-1);
        float lowerValue = cleanHood[lowerIndex];
        float upperValue = cleanHood[upperIndex];
        if(lowerIndex == upperIndex) {
            value = lowerValue;
        }
        else {
            assert(upperQuantile > lowerQuantile);
            assert(quantile >= lowerQuantile);
            float f = (quantile - lowerQuantile)/(upperQuantile - lowerQuantile);
            assert(f >= 0);
            assert(f <= 1);
            value   = lowerValue + (upperValue - lowerValue) * f;
        }
    }
    return value;
}
vec gridpp::util::calc_statistic(const vec2& array, gridpp::Statistic statistic) {
    int N = array.size();
    vec output(N);
    for(int n = 0; n < N; n++) {
        output[n] = gridpp::util::calc_statistic(array[n], statistic);
    }
    return output;
}
vec gridpp::util::calc_quantile(const vec2& array, float quantile) {
    int N = array.size();
    vec output(N);
    for(int n = 0; n < N; n++) {
        output[n] = gridpp::util::calc_quantile(array[n], quantile);
    }
    return output;
}
int gridpp::util::num_missing_values(const vec2& iArray) {
    int count = 0;
    for(int y = 0; y < iArray.size(); y++) {
        for(int x = 0; x < iArray[y].size(); x++) {
            count += !gridpp::util::is_valid(iArray[y][x]);
        }
    }
    return count;
}
void gridpp::util::debug(std::string string) {
    std::cout << string << std::endl;
}

void gridpp::util::warning(std::string string) {
    std::cout << "Warning: " << string << std::endl;
}

void gridpp::util::error(std::string iMessage) {
#ifdef DEBUG
    std::cout << "Error: " << iMessage << std::endl;
    void *array[10];
    size_t size = backtrace(array, 10);
    std::cout << "Stack trace:" << std::endl;
    backtrace_symbols_fd(array, size, 2);
    __gcov_flush();
#else
    std::cout << "Error: " << iMessage << std::endl;
#endif
    abort();
}
double gridpp::util::clock() {
    timeval t;
    gettimeofday(&t, NULL);
    double sec = (t.tv_sec);
    double msec= (t.tv_usec);
    return sec + msec/1e6;
}
vec gridpp::util::calc_even_quantiles(const vec& values, int num) {
    if(num >= values.size())
        return  values;

    vec quantiles;
    if(num == 0)
        return quantiles;

    quantiles.reserve(num);
    quantiles.push_back(values[0]);
    int count_lower = 0;
    for(int i = 0; i < values.size(); i++) {
        if(values[i] != quantiles[0])
            break;
        count_lower++;
    }
    int size = values.size();
    if(count_lower > size / num)
        quantiles.push_back(values[count_lower]);

    // Remove duplicates
    vec values_unique;
    values_unique.reserve(values.size());
    float last_quantile = quantiles[quantiles.size() - 1];
    for(int i = 0; i < values.size(); i++) {
        if(values[i] > last_quantile && (values_unique.size() == 0 || values[i] != values_unique[values_unique.size() - 1]))
            values_unique.push_back(values[i]);
    }
    if(values_unique.size() > 0) {
        int num_left = num - quantiles.size();
        for(int i = 1; i <= num_left; i++) {
            float f = float(i) / (num_left);
            int index = values_unique.size() * f - 1;
            if(index > 0) {
                float value = values_unique[index];
                quantiles.push_back(value);
            }
            else {
                std::cout << i << " " << f << " " << index << " " << num_left << " " << values_unique.size() << std::endl;
                std::cout << count_lower << " " << values.size() << " " << last_quantile << std::endl;
                throw std::runtime_error("Internal error in calc_even_quantiles.");
            }
        }
    }
    return quantiles;

}
int gridpp::util::get_lower_index(float iX, const std::vector<float>& iValues) {
    int index = gridpp::MV;
    for(int i = 0; i < (int) iValues.size(); i++) {
        float currValue = iValues[i];
        if(gridpp::util::is_valid(currValue)) {
            if(currValue < iX) {
                index = i;
            }
            else if(currValue == iX) {
                index = i;
                break;
            }
            else if(currValue > iX) {
                break;
            }
        }
    }
    return index;
}
int gridpp::util::get_upper_index(float iX, const std::vector<float>& iValues) {
    int index = gridpp::MV;
    for(int i = iValues.size()-1; i >= 0; i--) {
        float currValue = iValues[i];
        if(gridpp::util::is_valid(currValue)) {
            if(currValue > iX) {
                index = i;
            }
            else if(currValue == iX) {
                index = i;
                break;
            }
            else if(currValue < iX) {
                break;
            }
        }
    }
    return index;
}
float gridpp::util::interpolate(float x, const std::vector<float>& iX, const std::vector<float>& iY) {
    if(iX.size() != iY.size())
        throw std::invalid_argument("Dimension mismatch. Cannot interpolate.");
    float y = gridpp::MV;
    if(iX.size() == 0)
        return gridpp::MV;

    if(x >= iX[iX.size()-1])
        return iY[iX.size()-1];
    if(x <= iX[0])
        return iY[0];

    int i0   = get_lower_index(x, iX);
    int i1   = get_upper_index(x, iX);
    float x0 = iX[i0];
    float x1 = iX[i1];
    float y0 = iY[i0];
    float y1 = iY[i1];

    if(x0 == x1)
        y = (y0+y1)/2;
    else {
        assert(x1 >= x0);
        y = y0 + (y1 - y0) * (x - x0)/(x1 - x0);
    }

    return y;
}
void gridpp::util::not_implemented_error() {
    gridpp::util::error("Not implemented");
}
bool gridpp::util::compatible_size(const Grid& grid, const vec2& v) {
    return v.size() == 0 || grid.size()[0] != v.size() || grid.size()[1] != v[0].size();
}
bool gridpp::util::compatible_size(const Grid& grid, const vec3& v) {
    return v.size() == 0 || v[0].size() == 0 || grid.size()[0] != v[0].size() || grid.size()[1] != v[0][0].size();
}
vec2 gridpp::util::init_vec2(int Y, int X, float value) {
    vec2 output(Y);
    for(int y = 0; y < Y; y++)
        output[y].resize(X, value);
    return output;
}
float* gridpp::util::test_array(float* v, int n) {
    int count = 0;
    for(int i = 0; i < n; i++)
        count++;
    return v;
 }
float gridpp::util::test_vec_input(const vec& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++)
        total += input[i];
    return total;
}
int gridpp::util::test_ivec_input(const ivec& input) {
    int total = 0;
    for(int i = 0; i < input.size(); i++)
        total += input[i];
    return total;
}
float gridpp::util::test_vec2_input(const vec2& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++) {
        for(int j = 0; j < input[i].size(); j++) {
            total += input[i][j];
        }
    }
    return total;
}
float gridpp::util::test_vec3_input(const vec3& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++) {
        for(int j = 0; j < input[i].size(); j++) {
            for(int k = 0; k < input[i][j].size(); k++) {
                total += input[i][j][k];
            }
        }
    }
    return total;
}
vec gridpp::util::test_vec_output() {
    vec output(1000, 0);
    return output;
}
vec2 gridpp::util::test_vec2_output() {
    vec2 output(10);
    for(int i = 0; i < 10; i++)
        output[i].resize(10, 0);
    return output;
}
vec3 gridpp::util::test_vec3_output() {
    vec3 output(10);
    for(int i = 0; i < 10; i++) {
        output[i].resize(10);
        for(int j = 0; j < 10; j++)
            output[i][j].resize(10, 0);
    }
    return output;
}

ivec gridpp::util::test_ivec_output() {
    ivec output(1000, 0);
    return output;
}
ivec2 gridpp::util::test_ivec2_output() {
    ivec2 output(10);
    for(int i = 0; i < 10; i++)
        output[i].resize(10, 0);
    return output;
}

#if 0
vec2 gridpp::util::calc_gradient(const vec2& values, const vec2& aux, int radius) {
    int Y = values.size();
    int X = values[0].size();
    vec2 ret(values.size());

    for(int y = 0; y < Y; y++) {
        values[y].resize(X, 0);
        for(int x = 0; x < X; x++) {
            if(y > radius && y < Y - radius - 1 && x > radius && x < X - radius - 1) {
                vec xxx;
                vec yyy;
                xxx.reserve((2*radius+1)*(2*radius+1));
                yyy.reserve((2*radius+1)*(2*radius+1));
                for(int yy = y - radius; yy < Y - radius - 1; y++) {
                    for(int xx = x - radius; xx < X - radius - 1; x++) {
                        if(util::is_valid(aux[y][x])) {
                            xxx.push_back(aux[y][x]);
                            yyy.push_back(values[y][x]);
                        }
                    }
                }
                ivec coeffs = regression(xxx, yyy);
                ret[y][x] = coeffs[1];
            }
        }
    }
}
ivec gridpp::util::regression(const vec& x, const vec& y) {
    ivec ret(2, 0;
            float meanXY  = 0; // elev*T
            float meanX   = 0; // elev
            float meanY   = 0; // T
            float meanXX  = 0; // elev*elev
            int   counter = 0;
            int N = x.size();
            for(int n = 0; n < N; n++) {
            meanXY += x[n] * y[n];
            meanX += x[n];
            meanY += y[n];
            meanXX += x[n] * x[n];
            counter++;
            }
            meanXY /= counter;
            meanX  /= counter;
            meanY  /= counter;
            meanXX /= counter;
            float gradient = 0;
            if(meanXX - meanX*meanX != 0) {
            gradient = (meanXY - meanX*meanY)/(meanXX - meanX*meanX);
            }
            // TODO
            ret[0] = 0;
            ret[1] = gradient;
            return ret;
}
#endif
