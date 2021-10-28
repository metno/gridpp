#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec2 gridpp::calc_gradient(const vec2& base, const vec2& values, GradientType gradientType,
    int halfwidth, int min_num, float min_range, float default_gradient){

    if(halfwidth <= 0)
        throw std::invalid_argument("Halwidth cannot be <= 0; must be positive integer");
    if(min_range < 0)
        throw std::invalid_argument("min_range must be >= 0");
    if(base.size() == 0)
        throw std::invalid_argument("base input has no size");
    if(!gridpp::compatible_size(base, values))
        throw std::invalid_argument("base is not the same size as values");

    int nY = base.size();
    int nX = base[0].size();

    vec2 output = gridpp::init_vec2(nY, nX, default_gradient);

    // Min Max Gradient 
    if(gradientType == MinMax){
        for(int y = 0; y  < base.size(); y++) {
            for(int x = 0; x < base[y].size(); x++){
                float current_max = gridpp::MV;
                float current_min = gridpp::MV;
                int I_maxBase_Y = 0;
                int I_maxBase_X = 0;
                int I_minBase_Y = 0;
                int I_minBase_X = 0;
                for(int yy = std::max(0, y - halfwidth); yy <= std::min(nY - 1, y + halfwidth); yy++){
                    for(int xx = std::max(0, x - halfwidth); xx <= std::min(nX - 1, x + halfwidth); xx++){
                        float current_base = base[yy][xx];
                        if(!gridpp::is_valid(current_base)){
                            continue;
                        }
                        if(!gridpp::is_valid(values[yy][xx]))
                            continue;

                        if(!gridpp::is_valid(current_max) || current_base > current_max){
                            current_max = current_base;
                            I_maxBase_Y = yy;
                            I_maxBase_X = xx;
                        }
                        if(!gridpp::is_valid(current_min) || current_base < current_min){
                            current_min = current_base;
                            I_minBase_Y = yy;
                            I_minBase_X = xx;
                        }
                    }
                }
                if(!gridpp::is_valid(current_max) || !gridpp::is_valid(current_min)){
                    output[y][x] = default_gradient;
                }

                else if(abs(current_max - current_min) <= min_range){
                    output[y][x] = default_gradient;
                }
                else{
                    float diffBase = current_max - current_min;
                    float diffValues = values[I_maxBase_Y][I_maxBase_X] - values[I_minBase_Y][I_minBase_X];
                    output[y][x] = diffValues / diffBase;
                }
            }
        }
    }
    else if(gradientType == LinearRegression){
        vec2 base_x_base = gridpp::init_vec2(nY, nX, gridpp::MV);
        vec2 base_x_values = gridpp::init_vec2(nY, nX, gridpp::MV);

        vec2 meanX = gridpp::neighbourhood(base, halfwidth, gridpp::Mean);
        vec2 meanY = gridpp::neighbourhood(values, halfwidth, gridpp::Mean);

        #pragma omp parallel for collapse(2)
        for(int y = 0; y < nY; y++){
            for(int x = 0; x < nX; x++){
                base_x_base[y][x] = pow(base[y][x], 2);
                base_x_values[y][x] = base[y][x] * values[y][x];
            }
        }

        vec2 meanXX = gridpp::neighbourhood(base_x_base, halfwidth, gridpp::Mean);
        vec2 meanXY = gridpp::neighbourhood(base_x_values, halfwidth, gridpp::Mean);

        #pragma omp parallel for collapse(2)
        for(int y = 0; y < nY; y++){
            for(int x = 0; x < nX; x++){
                if(meanXX[y][x] - meanX[y][x] * meanX[y][x] != 0){
                    output[y][x] = (meanXY[y][x] - meanX[y][x] * meanY[y][x]) /
                                   (meanXX[y][x] - meanX[y][x] * meanX[y][x]);
                }
                else{
                    output[y][x] = default_gradient;
                }
            }
        }
    }
    return output;
}
