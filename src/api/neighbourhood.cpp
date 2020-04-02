#include "gridpp.h"

namespace {
    vec2 neighbourhood_brute_force(const vec2& input, int radius, std::string operation, float quantile);
    vec2 neighbourhood_quantile_ens(const vec3& input, float quantile, int radius, const vec& thresholds);
}


vec2 gridpp::neighbourhood_ens(const vec3& input, int iRadius, std::string operation) {
    vec2 flat(input.size());
    int Y = input.size();
    int X = input[0].size();
    int E = input[0][0].size();
    gridpp::util::StatType stat_type = gridpp::util::getStatType(operation);
    for(int y = 0; y < Y; y++) {
        flat[y].resize(X, 0);
    }
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            flat[y][x] = gridpp::util::calculate_stat(input[y][x], stat_type);
        }
    }
    vec2 results = neighbourhood(flat, iRadius, operation);
    return results;
}
vec2 gridpp::neighbourhood(const vec2& input, int iRadius, std::string operation) {
    double s_time = gridpp::util::clock();
    bool fast = true;
    int count_stat = 0;
    int nY = input.size();
    int nX = input[0].size();
    vec2 output(nY);
    gridpp::util::StatType stat_type = gridpp::util::getStatType(operation);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nX, 0);
    }
    if(stat_type == gridpp::util::StatTypeMean || stat_type == gridpp::util::StatTypeSum) {
        vec2 values;
        vec2 counts;
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
// #pragma omp parallel for
        for(int i = 0; i < nY; i++) {
            for(int j = 0; j < nX; j++) {
                int i1 = std::min(nY-1, i + iRadius);
                int j1 = std::min(nX-1, j + iRadius);
                int i0 = i - iRadius - 1;
                int j0 = j - iRadius - 1;
                float value11 = values[i1][j1];
                float value00 = 0;
                float value10 = 0;
                float value01 = 0;
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
                float value = value11 + value00 - value10 - value01;
                int count = count11 + count00 - count10 - count01;
                if(count > 0) {
                    if(stat_type == gridpp::util::StatTypeMean) {
                        value /= count;
                    }
                    output[i][j] = value;
                }
            }
        }
    }
    else if(gridpp::util::num_missing_values(input) == 0 && fast && (stat_type == gridpp::util::StatTypeMin || stat_type == gridpp::util::StatTypeMax)) {
        // Compute min/max quickly or any other quantile in a faster, but approximate way
        vec2 values;
        values.resize(nY);
        for(int i = 0; i < nY; i++) {
            values[i].resize(nX, 0);
        }
// #pragma omp parallel for
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
                    values[i][j] = gridpp::util::calculate_stat(neighbourhood, stat_type, 0);
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
                    slivers[j] = gridpp::util::calculate_stat(sliver, stat_type, 0);
                    count_stat += sliver.size();
                }
                for(int j = 0; j < nX; j++) {
                    std::vector<float> curr;
                    curr.reserve(2*iRadius);
                    for(int jj = std::max(0, j - iRadius); jj <= std::min(nX-1, j + iRadius); jj++) {
                        curr.push_back(slivers[jj]);
                    }
                    values[i][j] = gridpp::util::calculate_stat(curr, stat_type, 0);
                    count_stat += curr.size();
                }
            }
        }
// #pragma omp parallel for
        for(int i = 0; i < nY; i++) {
            for(int j = 0; j < nX; j++) {
                output[i][j] = values[i][j];
            }
        }
    }
    else {
        output = ::neighbourhood_brute_force(input, iRadius, operation, 0);
    }
    double e_time = gridpp::util::clock() ;
    // std::cout << count_stat << " " << e_time - s_time << " s" << std::endl;
    return output;
}
vec2 gridpp::neighbourhood_quantile(const vec2& input, float quantile, int radius, int num_thresholds) {
    vec3 input3(input.size());
    for(int i = 0; i < input.size(); i++) {
        input3[i].resize(input[i].size());
        for(int j = 0; j < input[i].size(); j++) {
            input3[i][j].push_back(input[i][j]);
        }
    }
    return gridpp::neighbourhood_quantile_ens(input3, quantile, radius, num_thresholds);
}
vec2 gridpp::neighbourhood_quantile_ens(const vec3& input, float quantile, int radius, int num_thresholds) {
    if(num_thresholds == 0) {
        // TODO
        abort();
    }
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

    return ::neighbourhood_quantile_ens(input, quantile, radius, thresholds);
}
namespace {
    vec2 neighbourhood_brute_force(const vec2& input, int iRadius, std::string operation, float quantile) {
        int count_stat = 0;
        int nY = input.size();
        int nX = input[0].size();
        vec2 output(nY);
        gridpp::util::StatType stat_type = gridpp::util::getStatType(operation);
        for(int y = 0; y < nY; y++) {
            output[y].resize(nX, 0);
        }
        vec2 values;
        values.resize(nY);
        for(int i = 0; i < nY; i++) {
            values[i].resize(nX, 0);
        }
// #pragma omp parallel for
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
                values[i][j] = gridpp::util::calculate_stat(neighbourhood, stat_type, quantile);
                count_stat += neighbourhood.size();
            }
        }
// #pragma omp parallel for
        for(int i = 0; i < nY; i++) {
            for(int j = 0; j < nX; j++) {
                output[i][j] = values[i][j];
            }
        }
    }
    vec2 neighbourhood_quantile_ens(const vec3& input, float quantile, int radius, const vec& thresholds) {
        double s_time = gridpp::util::clock();
        bool fast = true;
        assert(quantile >= 0);
        assert(quantile <= 1);
        int count_stat = 0;
        int nY = input.size();
        int nX = input[0].size();
        int nE = input[0][0].size();
        vec2 output(nY);

        for(int y = 0; y < nY; y++) {
            output[y].resize(nX, 0);
        }


        // Compute
        vec3 stats(thresholds.size());
        for(int t = 0; t < thresholds.size(); t++) {
            stats[t].resize(nY);
            for(int y = 0; y < nY; y++) {
                stats[t][y].resize(nX, 0);
           }
           vec2 temp(nY);
           for(int y = 0; y < nY; y++) {
               temp[y].resize(nX, 0);
               for(int x = 0; x < nX; x++) {
                   for(int e = 0; e < nE; e++) {
                       if(input[y][x][e] < thresholds[t]) {
                           temp[y][x] += 1.0 / nE;
                       }
                   }
               }
           }
           stats[t] = gridpp::neighbourhood(temp, radius, "mean");
        }
        for(int y = 0; y < nY; y++) {
            for(int x = 0; x < nX; x++) {
                vec yarray(thresholds.size());
                for(int t = 0; t < thresholds.size(); t++) {
                    for(int e = 0; e < nE; e++) {
                        yarray[t] += stats[t][y][x] / nE;
                    }
                }
                output[y][x] = gridpp::util::interpolate(quantile, yarray, thresholds);
                /*
                if(yarray[1] > 0.1) {
                    std::cout << output[y][x] << std::endl;
                    for(int t = 0; t < thresholds.size(); t++) {
                        std::cout << thresholds[t] << " " << stats[t][y][x] << std::endl;
                    }
                    abort();
                }
                */
            }
        }

        double e_time = gridpp::util::clock() ;
        std::cout << e_time - s_time << " s" << std::endl;
        return output;
    }
}
