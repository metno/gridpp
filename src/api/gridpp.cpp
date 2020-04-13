#include "gridpp.h"
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
