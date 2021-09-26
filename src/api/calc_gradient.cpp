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
        throw std::invalid_argument("Base input has no size");

    int nY = base.size();
    int nX = base[0].size();

    vec2 output = gridpp::init_vec2(base.size(), base[0].size(), default_gradient);

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
    // Linear Regression Gradient
    else if(gradientType == LinearRegression){
        for(int y = 0; y  < base.size(); y++) {
            for(int x = 0; x < base[y].size(); x++){
                float current_max = gridpp::MV;
                float current_min = gridpp::MV;

                float meanX = 0;
                float meanY = 0;
                float meanXY = 0;
                float meanXX = 0;

                int counter = 0;

                for(int yy = std::max(0, y - halfwidth); yy <= std::min(nY - 1, y + halfwidth); yy++){
                    for(int xx = std::max(0, x - halfwidth); xx <= std::min(nX - 1, x + halfwidth); xx++){
                        float current_base = base[yy][xx];
                        float current_values = values[yy][xx];
                        if(!gridpp::is_valid(current_max))
                            continue;
                        if(!gridpp::is_valid(current_values))
                            continue;

                        if(!gridpp::is_valid(current_max) || current_base > current_max){
                            current_max = current_base;
                        }

                        else if(!gridpp::is_valid(current_min) || current_base < current_min){
                            current_min = current_base;
                        }

                        meanX = meanX + current_base;
                        meanY = meanY + current_values;
                        meanXX = meanXX + current_base*current_base;
                        meanXY = meanXY + current_base*current_values;
                        counter = counter + 1;
                    }
                }

                if(counter == 0 || abs(current_max - current_min) <= min_range){
                    output[y][x] = default_gradient;
                }
                else {
                    meanX = meanX / counter;
                    meanY = meanY / counter;
                    meanXX = meanXX / counter;
                    meanXY = meanXY / counter;
                    if(meanXX - meanX*meanX != 0){
                        output[y][x] = (meanXY - meanX*meanY)/(meanXX - meanX*meanX);
                    }
                    else{
                        output[y][x] = default_gradient;
                    }
                }
            }
        }
    }
    return output;
}
