#include "gridpp.h"
#include <exception>
#include <iostream>

namespace {
    vec2 neighbourhood_brute_force(const vec2& input, int iRadius, gridpp::Statistic statistic, float quantile);
    vec2 neighbourhood_brute_force_ens(const vec3& input, int iRadius, gridpp::Statistic statistic, float quantile);
    vec3 vec2_to_vec3(const vec2& input);
}
vec2 gridpp::neighbourhood_ens(const vec3& input, int iRadius, gridpp::Statistic statistic) {
    vec2 flat(input.size());
    int Y = input.size();
    int X = input[0].size();
    int E = input[0][0].size();
    for(int y = 0; y < Y; y++) {
        flat[y].resize(X, 0);
    }
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            flat[y][x] = gridpp::util::calc_statistic(input[y][x], statistic);
        }
    }
    vec2 results = neighbourhood(flat, iRadius, statistic);
    return results;
}
vec2 gridpp::neighbourhood(const vec2& input, int iRadius, gridpp::Statistic statistic) {
    if(iRadius < 0)
        throw std::invalid_argument("Radius must be > 0");
    if(statistic == gridpp::Quantile)
        throw std::invalid_argument("Use neighbourhood_quantile for computing neighbourhood quantiles");
    if(input.size() == 0 || input[0].size() == 0)
        return vec2();

    double s_time = gridpp::util::clock();
    bool fast = true;
    int count_stat = 0;
    int nY = input.size();
    int nX = input[0].size();
    vec2 output(nY);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nX, gridpp::MV);
    }
    if(statistic == gridpp::Mean || statistic == gridpp::Sum) {
        dvec2 values;
        ivec2 counts;
        values.resize(nY);
        counts.resize(nY);
        for(int i = 0; i < nY; i++) {
            values[i].resize(nX, 0);
            counts[i].resize(nX, 0);
        }
        // Compute accumulated values
        for(int i = 0; i < nY; i++) {
            for(int j = 0; j < nX; j++) {
                float value = input[i][j];
                if(j == 0 && i == 0) {
                    // Lower corner
                    if(gridpp::util::is_valid(value)) {
                        values[i][j] = input[i][j];
                        counts[i][j] = 1;
                    }
                }
                else if(j == 0) {
                    // Lower row
                    if(gridpp::util::is_valid(value)) {
                        values[i][j] = values[i-1][j] + input[i][j];
                        counts[i][j] = counts[i-1][j] + 1;
                    }
                    else {
                        values[i][j] = values[i-1][j];
                        counts[i][j] = counts[i-1][j];
                    }
                }
                else if(i == 0) {
                    // Left column
                    if(gridpp::util::is_valid(value)) {
                        values[i][j] = values[i][j-1] + input[i][j];
                        counts[i][j] = counts[i][j-1] + 1;
                    }
                    else {
                        values[i][j] = values[i][j-1];
                        counts[i][j] = counts[i][j-1];
                    }

                }
                else {
                    if(gridpp::util::is_valid(value)) {
                        values[i][j] = values[i][j-1] + values[i-1][j] - values[i-1][j-1] + input[i][j];
                        counts[i][j] = counts[i][j-1] + counts[i-1][j] - counts[i-1][j-1] + 1;
                    }
                    else {
                        values[i][j] = values[i][j-1] + values[i-1][j] - values[i-1][j-1];
                        counts[i][j] = counts[i][j-1] + counts[i-1][j] - counts[i-1][j-1];
                    }
                }
            }
        }
        // Put neighbourhood into vector
        #pragma omp parallel for
        for(int i = 0; i < nY; i++) {
            for(int j = 0; j < nX; j++) {
                int i1 = std::min(nY-1, i + iRadius);
                int j1 = std::min(nX-1, j + iRadius);
                int i0 = i - iRadius - 1;
                int j0 = j - iRadius - 1;
                double value11 = values[i1][j1];
                double value00 = 0;
                double value10 = 0;
                double value01 = 0;
                int count11 = counts[i1][j1];
                int count00 = 0;
                int count10 = 0;
                int count01 = 0;
                if(i0 >= 0 && j0 >= 0) {
                    value00 = values[i0][j0];
                    value10 = values[i1][j0];
                    value01 = values[i0][j1];
                    count00 = counts[i0][j0];
                    count10 = counts[i1][j0];
                    count01 = counts[i0][j1];
                }
                else if(j0 >= 0) {
                    value10 = values[i1][j0];
                    count10 = counts[i1][j0];
                }
                else if(i0 >= 0) {
                    value01 = values[i0][j1];
                    count01 = counts[i0][j1];
                }
                double value = value11 + value00 - value10 - value01;
                int count = count11 + count00 - count10 - count01;
                if(count > 0) {
                    if(statistic == gridpp::Mean) {
                        value /= count;
                    }
                    output[i][j] = value;
                }
            }
        }
    }
    else if(statistic == gridpp::Min || statistic == gridpp::Max) {
        // Compute min/max quickly or any other quantile in a faster, but approximate way
        vec2 values;
        values.resize(nY);
        for(int i = 0; i < nY; i++) {
            values[i].resize(nX, 0);
        }
        #pragma omp parallel for
        for(int i = 0; i < nY; i++) {
            if(i < iRadius || i >= nY - iRadius) {
                // Regular way
                for(int j = 0; j < nX; j++) {
                    // Put neighbourhood into vector
                    std::vector<float> neighbourhood;
                    int Ni = std::min(nY-1, i+iRadius) - std::max(0, i-iRadius) + 1;
                    int Nj = std::min(nX-1, j+iRadius) - std::max(0, j-iRadius) + 1;
                    assert(Ni > 0);
                    assert(Nj > 0);
                    neighbourhood.reserve(Ni*Nj);
                    int index = 0;
                    for(int ii = std::max(0, i-iRadius); ii <= std::min(nY-1, i+iRadius); ii++) {
                        for(int jj = std::max(0, j-iRadius); jj <= std::min(nX-1, j+iRadius); jj++) {
                            float value = input[ii][jj];
                            if(gridpp::util::is_valid(value)) {
                                assert(index < Ni*Nj);
                                neighbourhood.push_back(value);
                                index++;
                            }
                        }
                    }
                    values[i][j] = gridpp::util::calc_statistic(neighbourhood, statistic);
                    count_stat += neighbourhood.size();
                }
            }
            else {
                // Fast way: Compute stats on each sliver
                std::vector<float> slivers(nX, 0);
                for(int j = 0; j < nX; j++) {
                    std::vector<float> sliver(2*iRadius+1, 0);
                    int count = 0;
                    for(int ii = i - iRadius; ii <= i + iRadius; ii++) {
                        sliver[count] = input[ii][j];
                        count++;
                    }
                    slivers[j] = gridpp::util::calc_statistic(sliver, statistic);
                    count_stat += sliver.size();
                }
                for(int j = 0; j < nX; j++) {
                    std::vector<float> curr;
                    curr.reserve(2*iRadius);
                    for(int jj = std::max(0, j - iRadius); jj <= std::min(nX-1, j + iRadius); jj++) {
                        curr.push_back(slivers[jj]);
                    }
                    values[i][j] = gridpp::util::calc_statistic(curr, statistic);
                    count_stat += curr.size();
                }
            }
        }
        #pragma omp parallel for
        for(int i = 0; i < nY; i++) {
            for(int j = 0; j < nX; j++) {
                output[i][j] = values[i][j];
            }
        }
    }
    else if(statistic == gridpp::Std || statistic == gridpp::Variance) {
        vec2 mean = gridpp::neighbourhood(input, iRadius, gridpp::Mean);
        vec2 input2(nY);
        for(int y = 0; y < nY; y++) {
            input2[y].resize(input[y].size(), 0);
            for(int x = 0; x < nX; x++) {
                input2[y][x] = input[y][x] * input[y][x];
            }
        }
        vec2 mean2 = gridpp::neighbourhood(input2, iRadius, gridpp::Mean);
        if(statistic == gridpp::Std) {
            for(int y = 0; y < nY; y++) {
                for(int x = 0; x < nX; x++) {
                    output[y][x] = sqrt(mean2[y][x] - mean[y][x] * mean[y][x]);
                }
            }
        }
        else {
            for(int y = 0; y < nY; y++) {
                for(int x = 0; x < nX; x++) {
                    output[y][x] = mean2[y][x] - mean[y][x] * mean[y][x];
                }
            }
        }
    }
    else {
        output = gridpp::neighbourhood_brute_force(input, iRadius, statistic);
    }
    double e_time = gridpp::util::clock() ;
    // std::cout << count_stat << " " << e_time - s_time << " s" << std::endl;
    return output;
}
vec gridpp::get_neighbourhood_thresholds(const vec2& input, int num_thresholds) {
    if(num_thresholds <= 0)
        throw std::invalid_argument("num_thresholds must be > 0");

    if(input.size() == 0 || input[0].size() == 0)
        return vec();

    int Y = input.size();
    assert(Y > 0);
    int X = input[0].size();
    assert(X > 0);
    size_t size = Y * X;
    vec all_values;
    all_values.reserve(size);
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            if(gridpp::util::is_valid(input[y][x])) {
                all_values.push_back(input[y][x]);
            }
        }
    }
    std::sort(all_values.begin(), all_values.end());
    vec thresholds = gridpp::util::calc_even_quantiles(all_values, num_thresholds);
    return thresholds;
}
vec gridpp::get_neighbourhood_thresholds(const vec3& input, int num_thresholds) {
    if(num_thresholds <= 0)
        throw std::invalid_argument("num_thresholds must be > 0");

    if(input.size() == 0 || input[0].size() == 0 || input[0][0].size() == 0)
        return vec();
    int Y = input.size();
    assert(Y > 0);
    int X = input[0].size();
    assert(X > 0);
    int E = input[0][0].size();
    assert(E > 0);
    size_t size = Y * X * E;
    vec all_values;
    all_values.reserve(size);
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            for(int e = 0; e < E; e++) {
                if(gridpp::util::is_valid(input[y][x][e])) {
                    all_values.push_back(input[y][x][e]);
                }
            }
        }
    }
    std::sort(all_values.begin(), all_values.end());
    vec thresholds = gridpp::util::calc_even_quantiles(all_values, num_thresholds);
    return thresholds;
}
vec2 gridpp::neighbourhood_quantile_fast(const vec2& input, float quantile, int radius, const vec& thresholds) {
    if(radius < 0)
        throw std::invalid_argument("Radius must be > 0");
    if(quantile < 0 || quantile > 1)
        throw std::invalid_argument("Quantile must be between 0 and 1 inclusive");

    if(input.size() == 0 || input[0].size() == 0)
        return vec2();
    double s_time = gridpp::util::clock();
    assert(quantile >= 0);
    assert(quantile <= 1);
    int nY = input.size();
    int nX = input[0].size();

    vec2 output(nY);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nX, gridpp::MV);
    }
    if(thresholds.size() == 0)
        return output;

    // Compute neighbourhood means for each threshold
    vec3 stats(thresholds.size());
    for(int t = 0; t < thresholds.size(); t++) {
        stats[t].resize(nY);
    }

    #pragma omp parallel for
    for(int t = 0; t < thresholds.size(); t++) {
        vec2 temp(nY);
        for(int y = 0; y < nY; y++) {
            temp[y].resize(nX, gridpp::MV);
            for(int x = 0; x < nX; x++) {
                int sum = 0;
                int count = 0;
                if(gridpp::util::is_valid(input[y][x])) {
                    if(input[y][x] <= thresholds[t]) {
                        sum++;
                    }
                    count++;
                }
                if(count > 0)
                    temp[y][x] = float(sum) / count;
            }
        }
        stats[t] = gridpp::neighbourhood(temp, radius, gridpp::Mean);
    }
    #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            vec yarray(thresholds.size(), gridpp::MV);
            bool is_missing = false;
            for(int t = 0; t < thresholds.size(); t++) {
                float sum = 0;
                int count = 0;
                if(gridpp::util::is_valid(stats[t][y][x])) {
                    sum = stats[t][y][x];
                    count++;
                }
                if(count > 0) {
                    yarray[t] = sum / count;
                    // Small floating point errors can occur in neighbourhood. Force values to be
                    // within [0, 1]
                    if(yarray[t] > 1)
                        yarray[t] = 1;
                    else if(yarray[t] < 0)
                        yarray[t] = 0;
                }
                else
                    is_missing = true;
            }
            if(!is_missing)
                output[y][x] = gridpp::util::interpolate(quantile, yarray, thresholds);
        }
    }

    double e_time = gridpp::util::clock() ;
    // std::cout << e_time - s_time << " s" << std::endl;
    return output;
}
vec2 gridpp::neighbourhood_quantile_ens_fast(const vec3& input, float quantile, int radius, const vec& thresholds) {
    if(radius < 0)
        throw std::invalid_argument("Radius must be > 0");
    if(quantile < 0 || quantile > 1)
        throw std::invalid_argument("Quantile must be between 0 and 1 inclusive");

    if(input.size() == 0 || input[0].size() == 0 || input[0][0].size() == 0)
        return vec2();
    double s_time = gridpp::util::clock();
    assert(quantile >= 0);
    assert(quantile <= 1);
    int nY = input.size();
    int nX = input[0].size();
    int nE = input[0][0].size();

    vec2 output(nY);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nX, gridpp::MV);
    }
    if(thresholds.size() == 0)
        return output;

    // Compute neighbourhood means for each threshold
    vec3 stats(thresholds.size());
    for(int t = 0; t < thresholds.size(); t++) {
        stats[t].resize(nY);
    }

    #pragma omp parallel for
    for(int t = 0; t < thresholds.size(); t++) {
        vec2 temp(nY);
        for(int y = 0; y < nY; y++) {
            temp[y].resize(nX, gridpp::MV);
            for(int x = 0; x < nX; x++) {
                int sum = 0;
                int count = 0;
                for(int e = 0; e < nE; e++) {
                    if(gridpp::util::is_valid(input[y][x][e])) {
                        if(input[y][x][e] <= thresholds[t]) {
                            sum++;
                        }
                        count++;
                    }
                }
                if(count > 0)
                    temp[y][x] = float(sum) / count;
            }
        }
        stats[t] = gridpp::neighbourhood(temp, radius, gridpp::Mean);
    }
    #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            vec yarray(thresholds.size(), gridpp::MV);
            bool is_missing = false;
            for(int t = 0; t < thresholds.size(); t++) {
                float sum = 0;
                int count = 0;
                for(int e = 0; e < nE; e++) {
                    if(gridpp::util::is_valid(stats[t][y][x])) {
                        sum += stats[t][y][x];
                        count++;
                    }
                }
                if(count > 0) {
                    yarray[t] = sum / count;
                    // Small floating point errors can occur in neighbourhood. Force values to be
                    // within [0, 1]
                    if(yarray[t] > 1)
                        yarray[t] = 1;
                    else if(yarray[t] < 0)
                        yarray[t] = 0;
                }
                else
                    is_missing = true;
            }
            if(!is_missing)
                output[y][x] = gridpp::util::interpolate(quantile, yarray, thresholds);
        }
    }

    double e_time = gridpp::util::clock() ;
    // std::cout << e_time - s_time << " s" << std::endl;
    return output;
}
vec2 gridpp::neighbourhood_brute_force(const vec2& input, int iRadius, gridpp::Statistic statistic) {
    return ::neighbourhood_brute_force(input, iRadius, statistic, 0);
}
vec2 gridpp::neighbourhood_quantile(const vec2& input, float quantile, int iRadius) {
    return ::neighbourhood_brute_force(input, iRadius, gridpp::Quantile, quantile);
}
vec2 gridpp::neighbourhood_quantile_ens(const vec3& input, float quantile, int iRadius) {
    return ::neighbourhood_brute_force_ens(input, iRadius, gridpp::Quantile, quantile);
}
namespace {
    vec2 neighbourhood_brute_force(const vec2& input, int iRadius, gridpp::Statistic statistic, float quantile) {
        if(iRadius < 0)
            throw std::invalid_argument("Radius must be > 0");
        if(input.size() == 0 || input[0].size() == 0)
            return vec2();

        int count_stat = 0;
        int nY = input.size();
        if(nY == 0)
            return vec2(0);
        int nX = input[0].size();
        if(nX == 0)
            return vec2(0);
        vec2 output(nY);
        for(int y = 0; y < nY; y++) {
            output[y].resize(nX, gridpp::MV);
        }
        #pragma omp parallel for
        for(int i = 0; i < nY; i++) {
            for(int j = 0; j < nX; j++) {
                // Put neighbourhood into vector
                std::vector<float> neighbourhood;
                int Ni = std::min(nY-1, i+iRadius) - std::max(0, i-iRadius) + 1;
                int Nj = std::min(nX-1, j+iRadius) - std::max(0, j-iRadius) + 1;
                assert(Ni > 0);
                assert(Nj > 0);
                neighbourhood.resize(Ni*Nj, gridpp::MV);
                int index = 0;
                for(int ii = std::max(0, i-iRadius); ii <= std::min(nY-1, i+iRadius); ii++) {
                    for(int jj = std::max(0, j-iRadius); jj <= std::min(nX-1, j+iRadius); jj++) {
                        float value = input[ii][jj];
                        assert(index < Ni*Nj);
                        neighbourhood[index] = value;
                        index++;
                    }
                }
                assert(index == Ni*Nj);
                if(statistic == gridpp::Quantile)
                    output[i][j] = gridpp::util::calc_quantile(neighbourhood, quantile);
                else
                    output[i][j] = gridpp::util::calc_statistic(neighbourhood, statistic);
                count_stat += neighbourhood.size();
            }
        }
        return output;
    }
    vec2 neighbourhood_brute_force_ens(const vec3& input, int iRadius, gridpp::Statistic statistic, float quantile) {
        if(iRadius < 0)
            throw std::invalid_argument("Radius must be > 0");
        if(input.size() == 0 || input[0].size() == 0 || input[0][0].size() == 0)
            return vec2();

        int count_stat = 0;
        int nY = input.size();
        if(nY == 0)
            return vec2();
        int nX = input[0].size();
        if(nX == 0)
            return vec2();
        int nEns = input[0][0].size();
        if(nEns == 0)
            return vec2();
        vec2 output(nY);
        for(int y = 0; y < nY; y++) {
            output[y].resize(nX, gridpp::MV);
        }
        #pragma omp parallel for
        for(int i = 0; i < nY; i++) {
            for(int j = 0; j < nX; j++) {
                // Put neighbourhood into vector
                std::vector<float> neighbourhood;
                int Ni = std::min(nY-1, i+iRadius) - std::max(0, i-iRadius) + 1;
                int Nj = std::min(nX-1, j+iRadius) - std::max(0, j-iRadius) + 1;
                assert(Ni > 0);
                assert(Nj > 0);
                neighbourhood.resize(Ni*Nj*nEns, gridpp::MV);
                int index = 0;
                for(int ii = std::max(0, i-iRadius); ii <= std::min(nY-1, i+iRadius); ii++) {
                    for(int jj = std::max(0, j-iRadius); jj <= std::min(nX-1, j+iRadius); jj++) {
                        for(int e = 0; e < nEns; e++) {
                            float value = input[ii][jj][e];
                            assert(index < Ni*Nj);
                            neighbourhood[index] = value;
                            index++;
                        }
                    }
                }
                assert(index == Ni*Nj);
                if(statistic == gridpp::Quantile)
                    output[i][j] = gridpp::util::calc_quantile(neighbourhood, quantile);
                else
                    output[i][j] = gridpp::util::calc_statistic(neighbourhood, statistic);
                count_stat += neighbourhood.size();
            }
        }
        return output;
    }
    vec3 vec2_to_vec3(const vec2& input) {
        if(input.size() == 0 || input[0].size() == 0)
            return vec3();
        int Y = input.size();
        int X = input[0].size();
        vec3 output(Y);
        for(int y = 0; y < Y; y++) {
            output[y].resize(X);
            for(int x = 0; x < X; x++) {
                output[y][x].push_back(input[y][x]);
            }
        }
        return output;
    }
}
