#include "gridpp.h"
#include <cstdlib>
#include <string>
#include <sstream>

using namespace gridpp;

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
    else if(name == "count") {
        type = gridpp::Count;
    }
    else if(name == "randomchoice") {
        type = gridpp::RandomChoice;
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
int gridpp::get_omp_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#endif
    return 0;
}

void gridpp::set_debug_level(int level) {
    gridpp::debug_level = level;
}

int gridpp::get_debug_level() {
    return gridpp::debug_level;
}
gridpp::not_implemented_exception::not_implemented_exception() :
    std::logic_error("Function not yet implemented") {
};
