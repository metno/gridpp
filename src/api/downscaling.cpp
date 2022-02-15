#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

vec gridpp::downscaling(const Grid& igrid, const Points& opoints, const vec2& ivalues, Downscaler downscaler) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");

    vec output;
    if(downscaler == gridpp::Nearest)
        output = gridpp::nearest(igrid, opoints, ivalues);
    else if(downscaler == gridpp::Bilinear)
        output = gridpp::bilinear(igrid, opoints, ivalues);
    else
        throw std::invalid_argument("Invalid downscaler");
    return output;
}

vec2 gridpp::downscaling(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, Downscaler downscaler) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");

    vec2 output;
    if(downscaler == gridpp::Nearest)
        output = gridpp::nearest(igrid, ogrid, ivalues);
    else if(downscaler == gridpp::Bilinear)
        output = gridpp::bilinear(igrid, ogrid, ivalues);
    else
        throw std::invalid_argument("Invalid downscaler");
    return output;
}

vec2 gridpp::downscaling(const Grid& igrid, const Points& opoints, const vec3& ivalues, Downscaler downscaler) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");

    vec2 output;
    if(downscaler == gridpp::Nearest)
        output = gridpp::nearest(igrid, opoints, ivalues);
    else if(downscaler == gridpp::Bilinear)
        output = gridpp::bilinear(igrid, opoints, ivalues);
    else
        throw std::invalid_argument("Invalid downscaler");
    return output;
}

vec3 gridpp::downscaling(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, Downscaler downscaler) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");

    vec3 output;
    if(downscaler == gridpp::Nearest)
        output = gridpp::nearest(igrid, ogrid, ivalues);
    else if(downscaler == gridpp::Bilinear)
        output = gridpp::bilinear(igrid, ogrid, ivalues);
    else
        throw std::invalid_argument("Invalid downscaler");
    return output;
}
