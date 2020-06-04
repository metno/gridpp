#include "gridpp.h"
#include <cstdlib>
#include <string>
#include <sstream>

std::string gridpp::version() {
    return __version__;
}
gridpp::Statistic gridpp::get_statistic(std::string name) {
    gridpp::Statistic type;
    if(name == "mean") {
        type = gridpp::Mean;
    }
    else if(name == "min") {
        type = gridpp::Min;
    }
    else if(name == "max") {
        type = gridpp::Max;
    }
    else if(name == "median") {
        type = gridpp::Median;
    }
    else if(name == "quantile") {
        type = gridpp::Quantile;
    }
    else if(name == "std") {
        type = gridpp::Std;
    }
    else if(name == "sum") {
        type = gridpp::Sum;
    }
    else {
        type = gridpp::Unknown;
    }
    return type;
}
void gridpp::initialize_omp() {
#ifdef _OPENMP
    int num_threads = 1;
    const char* num_threads_char = std::getenv("OMP_NUM_THREADS");
    if(num_threads_char != NULL) {
        std::istringstream(std::string(num_threads_char)) >> num_threads;
        if(num_threads <= 0)
            num_threads = 1;
    }
    gridpp::set_omp_threads(num_threads);
#endif
}
void gridpp::set_omp_threads(int num) {
#ifdef _OPENMP
    // omp_set_dynamic(0);
    omp_set_num_threads(num);
#endif
}

float* gridpp::test_array(float* v, int n) {
    int count = 0;
    for(int i = 0; i < n; i++)
        count++;
    return v;
 }
int gridpp::test_vec_input(const vec& input) {
    return 0;
}
int gridpp::test_vec2_input(const vec2& input) {
    return 0;
}
int gridpp::test_vec3_input(const vec3& input) {
    return 0;
}
vec1 gridpp::test_vec_output() {
    vec output(1000, 0);
    return output;
}
vec2 gridpp::test_vec2_output() {
    vec2 output(1000);
    for(int i = 0; i < 1000; i++)
        output[i].resize(1000, 0);
    return output;
}
vec3 gridpp::test_vec3_output() {
    vec3 output(1000);
    for(int i = 0; i < 1000; i++) {
        output[i].resize(1000);
        for(int j = 0; j < 1000; j++)
            output[i][j].resize(1000, 0);
    }
    return output;
}
