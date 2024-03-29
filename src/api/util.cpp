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

using namespace gridpp;

bool gridpp::is_valid(float value) {
    return !std::isnan(value) && !std::isinf(value) && value != gridpp::MV;
}
float gridpp::calc_statistic(const vec& array, gridpp::Statistic statistic) {
    // Initialize to missing
    float value = gridpp::MV;
    if(statistic == gridpp::Mean || statistic == gridpp::Sum || statistic == gridpp::Count) {
        float total = 0;
        int count = 0;
        for(int n = 0; n < array.size(); n++) {
            if(gridpp::is_valid(array[n])) {
                total += array[n];
                count++;
            }
        }
        if (statistic == gridpp::Count)
            value = count;
        else if(count > 0) {
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
            if(gridpp::is_valid(array[n])) {
                if(!gridpp::is_valid(K))
                    K = array[n];
                assert(gridpp::is_valid(K));

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
    else if (statistic == gridpp::RandomChoice) {
        int num_valid = 0;
        for(int n = 0; n < array.size(); n++) {
            if(gridpp::is_valid(array[n])) {
                num_valid++;
            }
        }
        if(num_valid > 0) {
            int random_index = rand() % num_valid;
            assert(random_index < array.size());
            int count = 0;
            for(int n = 0; n < array.size(); n++) {
                if(gridpp::is_valid(array[n])) {
                    if(count == random_index) {
                        value = array[n];
                        break;
                    }
                    count++;
                }
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
        value = gridpp::calc_quantile(array, quantile);
    }
    return value;
}
float gridpp::calc_quantile(const vec& array, float quantile) {
    int T = array.size();
    if(quantile < 0 || quantile > 1) {
        throw std::invalid_argument("calc_quantile: Quantile must be between 0 and 1 inclusive");
    }
    if(!gridpp::is_valid(quantile))
        return gridpp::MV;

    if(T == 0)
        return gridpp::MV;
    if(quantile == 0) {
        float min = gridpp::MV;
        for(int i = 0; i < T; i++) {
            float val = array[i];
            if(!gridpp::is_valid(val))
                continue;
            else if(!gridpp::is_valid(min))
                min = val;
            else if(val < min)
                min = val;
        }
        return min;
    }
    else if(quantile == 1) {
        float max = gridpp::MV;
        for(int i = 0; i < T; i++) {
            float val = array[i];
            if(!gridpp::is_valid(val))
                continue;
            else if(!gridpp::is_valid(max))
                max = val;
            else if(val > max)
                max = val;
        }
        return max;
    }
    // Initialize to missing
    float value = gridpp::MV;
    // Remove missing
    vec cleanHood;
    cleanHood.reserve(array.size());
    for(int i = 0; i < array.size(); i++) {
        if(gridpp::is_valid(array[i]))
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
vec gridpp::calc_quantile(const vec2& array, float quantile) {
    int N = array.size();
    vec output(N);
    for(int n = 0; n < N; n++) {
        output[n] = gridpp::calc_quantile(array[n], quantile);
    }
    return output;
}
vec2 gridpp::calc_quantile(const vec3& array, const vec2& quantile) {
    if(!gridpp::compatible_size(quantile, array))
        throw std::invalid_argument("Dimension mismatch between array and quantile");
    int Y = array.size();
    if(Y == 0)
        return vec2();
    int X = array[0].size();
    if(X == 0)
        return vec2();
    int T = array[0][0].size();
    if(T == 0) {
        vec2 output = gridpp::init_vec2(Y, X, gridpp::MV);
        return output;
    }
    vec2 output = gridpp::init_vec2(Y, X);
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            output[y][x] = gridpp::calc_quantile(array[y][x], quantile[y][x]);
        }
    }
    return output;
}
vec gridpp::calc_statistic(const vec2& array, gridpp::Statistic statistic) {
    int N = array.size();
    vec output(N);
    for(int n = 0; n < N; n++) {
        output[n] = gridpp::calc_statistic(array[n], statistic);
    }
    return output;
}
int gridpp::num_missing_values(const vec2& iArray) {
    int count = 0;
    for(int y = 0; y < iArray.size(); y++) {
        for(int x = 0; x < iArray[y].size(); x++) {
            count += !gridpp::is_valid(iArray[y][x]);
        }
    }
    return count;
}
void gridpp::debug(std::string string) {
    std::cout << string << std::endl;
}

void gridpp::warning(std::string string) {
    std::cout << "Warning: " << string << std::endl;
}

void gridpp::error(std::string iMessage) {
#ifdef DEBUG
    std::cout << "Error: " << iMessage << std::endl;
    void *array[10];
    size_t size = backtrace(array, 10);
    std::cout << "Stack trace:" << std::endl;
    backtrace_symbols_fd(array, size, 2);
#else
    std::cout << "Error: " << iMessage << std::endl;
#endif
    throw std::runtime_error(iMessage);
}
void gridpp::future_deprecation_warning(std::string function, std::string other) {
    std::cout << "Future deprecation warning: " << function << " will be deprecated";
    if(other != "")
        std::cout << ", use " << other << " instead." << std::endl;
    else
        std::cout << "." << std::endl;
}

double gridpp::clock() {
    timeval t;
    gettimeofday(&t, NULL);
    double sec = (t.tv_sec);
    double msec= (t.tv_usec);
    return sec + msec/1e6;
}
vec gridpp::calc_even_quantiles(const vec& values, int num) {
    vec quantiles;
    int size = values.size();
    if(num == 0 || size == 0)
        return quantiles;

    vec sorted = values;
    std::sort(sorted.begin(), sorted.end());

    // Asking for more quantiles than we have data points
    if(num >= size) {
        // Find the unique values
        quantiles.reserve(num);
        quantiles.push_back(sorted[0]);
        for(int i = 1; i < size; i++) {
            if(sorted[i] != sorted[i-1])
                quantiles.push_back(sorted[i]);
        }
        return  quantiles;
    }


    float lowest = sorted[0];
    float highest = sorted[size - 1];

    // If there are multiple identical values at the start, pick the first value past the set of
    // identical values
    int count_lower = 0;
    for(int i = 0; i < sorted.size(); i++) {
        if(sorted[i] != lowest)
            break;
        count_lower++;
    }
    quantiles.reserve(num);
    quantiles.push_back(lowest);
    float last_added = lowest;
    if(num == 2) {
        if(lowest != highest)
            quantiles.push_back(highest);
        return quantiles;
    }

    // Add the first past the lowest set of repeated values
    bool repeated_at_beginning = count_lower < size && count_lower > size / num;
    if(repeated_at_beginning) {
        assert(count_lower < sorted.size());
        quantiles.push_back(sorted[count_lower]);
    }
    last_added = quantiles[quantiles.size() - 1];

    // Remove duplicates
    vec remaining_unique_values;
    remaining_unique_values.reserve(sorted.size());
    for(int i = 0; i < sorted.size(); i++) {
        if(sorted[i] > last_added && (remaining_unique_values.size() == 0 || sorted[i] != remaining_unique_values[remaining_unique_values.size() - 1]))
            remaining_unique_values.push_back(sorted[i]);
    }

    if(remaining_unique_values.size() > 0) {
        int num_left = num - quantiles.size();
        for(int i = 1; i <= num_left; i++) {
            float f = float(i) / (num_left);
            int index = remaining_unique_values.size() * f - 1;
            if(index >= 0) {
                assert(index < remaining_unique_values.size());
                float value = remaining_unique_values[index];
                quantiles.push_back(value);
            }
            else {
                std::cout << i << " " << f << " " << index << " " << num_left << " " << remaining_unique_values.size() << std::endl;
                std::cout << count_lower << " " << sorted.size() << " " << last_added << std::endl;
                throw std::runtime_error("Internal error in calc_even_quantiles.");
            }
        }
    }
    return quantiles;

}
int gridpp::get_lower_index(float iX, const vec& iValues) {
    int index = gridpp::MV;
    for(int i = 0; i < (int) iValues.size(); i++) {
        float currValue = iValues[i];
        if(gridpp::is_valid(currValue)) {
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
int gridpp::get_upper_index(float iX, const vec& iValues) {
    int index = gridpp::MV;
    for(int i = iValues.size()-1; i >= 0; i--) {
        float currValue = iValues[i];
        if(gridpp::is_valid(currValue)) {
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
float gridpp::interpolate(float x, const vec& iX, const vec& iY) {
    if(!gridpp::is_valid(x))
        return gridpp::MV;
    if(iX.size() != iY.size())
        throw std::invalid_argument("Dimension mismatch. Cannot interpolate.");
    float y = gridpp::MV;
    if(iX.size() == 0)
        return gridpp::MV;

    if(x > iX[iX.size()-1])
        return iY[iX.size()-1];
    if(x < iX[0])
        return iY[0];

    int i0   = get_lower_index(x, iX);
    int i1   = get_upper_index(x, iX);
    float x0 = iX[i0];
    float x1 = iX[i1];
    float y0 = iY[i0];
    float y1 = iY[i1];

    if(x0 == x1) {
        if(i0 == 0 && i1 == iX.size()-1)
            y = (y0 + y1)/2;
        else if(i0 == 0)
            y = y1;
        else if(i1 == iX.size()-1)
            y = y0;
        else
            y = (y0+y1)/2;
    }
    else {
        assert(x1 >= x0);
        y = y0 + (y1 - y0) * (x - x0)/(x1 - x0);
    }

    return y;
}
vec gridpp::interpolate(const vec& x, const vec& iX, const vec& iY) {
    if(iX.size() != iY.size())
        throw std::invalid_argument("Dimension mismatch. Cannot interpolate.");

    vec y(x.size());

    #pragma omp parallel for
    for(int i = 0; i < x.size(); i++) {
        y[i] = interpolate(x[i], iX, iY);
    }
    return y;
}
bool gridpp::compatible_size(const Grid& grid, const vec2& v) {
    return v.size() == 0 || (grid.size()[0] == v.size() && grid.size()[1] == v[0].size());
}
bool gridpp::compatible_size(const Grid& grid, const vec3& v) {
    return v.size() == 0 || (v[0].size() == 0 || (grid.size()[0] == v[0].size() && grid.size()[1] == v[0][0].size()));
}
bool gridpp::compatible_size(const Points& points, const vec2& v) {
    return v.size() == 0 || points.size() == v[0].size();
}
bool gridpp::compatible_size(const Points& points, const vec& v) {
    return points.size() == v.size();
}
bool gridpp::compatible_size(const vec2& a, const vec2& b) {
    if(a.size() != b.size())
        return false;
    for(int y = 0; y < a.size(); y++) {
        if(a[y].size() != b[y].size())
            return false;
    }
    return true;
}
bool gridpp::compatible_size(const vec2& a, const vec3& b) {
    // No y
    if(a.size() == 0 && b.size() == 0)
        return true;
    if(a.size() != b.size())
        return false;

    // No x
    if(a[0].size() == 0 && b[0].size() == 0)
        return true;
    if(a[0].size() != b[0].size())
        return false;
    return true;
}
bool gridpp::compatible_size(const vec3& a, const vec3& b) {
    if(a.size() != b.size())
        return false;
    for(int t = 0; t < a.size(); t++) {
        if(a[t].size() != b[t].size())
            return false;
        for(int y = 0; y < a[t].size(); y++) {
            if(a[t][y].size() != b[t][y].size())
                return false;
        }
    }
    return true;
}
vec2 gridpp::init_vec2(int Y, int X, float value) {
    vec2 output(Y);
    for(int y = 0; y < Y; y++)
        output[y].resize(X, value);
    return output;
}
ivec2 gridpp::init_ivec2(int Y, int X, int value) {
    ivec2 output(Y);
    for(int y = 0; y < Y; y++)
        output[y].resize(X, value);
    return output;
}
vec3 gridpp::init_vec3(int Y, int X, int T, float value) {
    vec3 output(Y);
    for(int y = 0; y < Y; y++) {
        output[y].resize(X);
        for(int x = 0; x < X; x++)
            output[y][x].resize(T, value);
    }
    return output;
}
ivec3 gridpp::init_ivec3(int Y, int X, int T, int value) {
    ivec3 output(Y);
    for(int y = 0; y < Y; y++) {
        output[y].resize(X);
        for(int x = 0; x < X; x++)
            output[y][x].resize(T, value);
    }
    return output;
}
#if 0
vec2 gridpp::calc_gradient(const vec2& values, const vec2& aux, int radius) {
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
ivec gridpp::regression(const vec& x, const vec& y) {
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
namespace {
    Point vect2d(const Point& p1, const Point& p2) {
        Point temp(0, 0);
        temp.lon = (p2.lon - p1.lon);
        temp.lat = -1 * (p2.lat - p1.lat);
        return temp;
    }
}
// Based on solution posted on Stack Overflow: https://stackoverflow.com/a/42396910/2540278
bool gridpp::point_in_rectangle(const Point& A, const Point& B, const Point& C, const Point& D, const Point& m) {
    Point AB = ::vect2d(A, B);  float C1 = -1 * (AB.lat*A.lon + AB.lon*A.lat); float  D1 = (AB.lat*m.lon + AB.lon*m.lat) + C1;
    Point AD = ::vect2d(A, D);  float C2 = -1 * (AD.lat*A.lon + AD.lon*A.lat); float D2 = (AD.lat*m.lon + AD.lon*m.lat) + C2;
    Point BC = ::vect2d(B, C);  float C3 = -1 * (BC.lat*B.lon + BC.lon*B.lat); float D3 = (BC.lat*m.lon + BC.lon*m.lat) + C3;
    Point CD = ::vect2d(C, D);  float C4 = -1 * (CD.lat*C.lon + CD.lon*C.lat); float D4 = (CD.lat*m.lon + CD.lon*m.lat) + C4;
    // std::cout << " # " << D1 << " " << D2 << " " << D3 << " " << D4 << std::endl;
    // Clockwise and counter-clockwise alternatives
    bool opt1 = 0 >= D1 && 0 >= D4 && 0 <= D2 && 0 >= D3;
    bool opt2 = 0 <= D1 && 0 <= D4 && 0 >= D2 && 0 <= D3;
    return opt1 || opt2;
}

bool gridpp::convert_coordinates(const vec& lats, const vec& lons, CoordinateType type, vec& x_coords, vec& y_coords, vec& z_coords) {
    int N = lats.size();
    x_coords.resize(N);
    y_coords.resize(N);
    z_coords.resize(N);
    for(int i = 0; i < N; i++) {
        convert_coordinates(lats[i], lons[i], type, x_coords[i], y_coords[i], z_coords[i]);
    }

    return true;
}

bool gridpp::convert_coordinates(float lat, float lon, CoordinateType type, float& x_coord, float& y_coord, float& z_coord) {
    if(!gridpp::is_valid_lat(lat, type) || !gridpp::is_valid_lon(lon, type)) {
        std::stringstream ss;
        ss << "Invalid coords: " << lat << "," << lon << std::endl;
        throw std::invalid_argument(ss.str());
    }
    if(type == gridpp::Cartesian) {
        x_coord = lon;
        y_coord = lat;
        z_coord = 0;
    }
    else {
        double lonr = M_PI / 180 * lon;
        double latr = M_PI / 180 * lat;
        // std::cout << lon << " " << lat << std::endl;
        x_coord = std::cos(latr) * std::cos(lonr) * gridpp::radius_earth;
        y_coord = std::cos(latr) * std::sin(lonr) * gridpp::radius_earth;
        z_coord = std::sin(latr) * gridpp::radius_earth;
    }
    return true;
}

bool gridpp::is_valid_lat(float lat, CoordinateType type) {
    if(type == gridpp::Cartesian)
        return gridpp::is_valid(lat);
    return gridpp::is_valid(lat) && (lat >= -90.001) && (lat <= 90.001);
};
bool gridpp::is_valid_lon(float lon, CoordinateType type) {
    return gridpp::is_valid(lon);
}
