#include "gridpp.h"
std::string gridpp::version() {
    return __version__;
}

gridpp::Operation gridpp::get_operation(std::string name) {
    gridpp::Operation type;
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
    return type;
}
