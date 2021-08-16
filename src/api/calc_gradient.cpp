#include "gridpp.h"
#include <iostream>

using namespace gridpp;



/*GradientType{
*   minMax = 0;
*   LinearRegression = 10;
};*/

vec2 gridpp::calc_gradient(const vec2& base, const vec2& values, //GradientType gradientType,
    int halfwidth , int min_num, float min_range, float default_gradient){

    if(halfwidth < 0)
        throw std::invalid_argument("halfwidth must be positive");

    std::cout << " Starting ";

    int gradientType = 0;

    //if(halfwidth == 0)
    //    throw std::invalid_argument("Halwidth cannot be 0; must be positive integer");
    //if(base.size() == 0)
    //    throw std::invalid_argument("Base input has no size");

    // calculate maximum element and minimum element within halfwidth of each [i][j]
    int nY = base.size();
    int nX = base[0].size();

    vec2 output = gridpp::init_vec2(base.size(), base[0].size());

    // Min Max Gradient 
    if(gradientType == 0){
        for(int y = 0; y  < base.size(); y++) {
            for(int x = 0; x < base[y].size(); x++){
                

                bool start = true;

                float current_max = 0;
                float current_min = 0;
                int I_maxBase_Y = 0;
                int I_maxBase_X = 0;
                int I_minBase_Y = 0;
                int I_minBase_X = 0;
                
                for(int yy = std::max(0, y - halfwidth); yy <= std::min(nY - 1, y + halfwidth); yy++){               
                    for(int xx = std::max(0, x - halfwidth); xx <= std::min(nX - 1, x + halfwidth); xx++){

                        float current_base = base[yy][xx];
                        if(start){
                            start = false;
                            current_max = current_base;
                            I_maxBase_Y = yy;
                            I_maxBase_X = xx;
                            current_min = current_base;
                            I_minBase_Y = yy;
                            I_minBase_X = xx;
                        }

                        else if(current_base > current_max){
                            current_max = current_base;
                            I_maxBase_Y = yy;
                            I_maxBase_X = xx;
                        }
                        else if(current_base < current_min){
                            current_min = current_base;
                            I_minBase_Y = yy;
                            I_minBase_X = xx;
                        }
                    }       
                }
                if(abs(current_max - current_min) < min_range){
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
    else if(gradientType == 10){
        for(int y = 0; y  < base.size(); y++) {
            for(int x = 0; x < base[y].size(); x++){             
                bool start = true;
                float current_max = 0;
                float current_min = 0;

                float meanX = 0; 
                float meanY = 0;
                float meanXY = 0;
                float meanXX = 0;

                int counter = 0;
                
                for(int yy = std::max(0, y - halfwidth); yy <= std::min(nY - 1, y + halfwidth); yy++){               
                    for(int xx = std::max(0, x - halfwidth); xx <= std::min(nX - 1, x + halfwidth); xx++){
                        float current_base = base[yy][xx];
                        float current_values = values[yy][xx];

                        if(start){
                            start = false;
                            current_max = current_base;
                            current_min = current_base;
                        }

                        else if(current_base > current_max){
                            current_max = current_base;
                        }

                        else if(current_base < current_min){
                            current_min = current_base;
                        }

                        meanX = meanX + current_base;
                        meanY = meanY + current_values;
                        meanXX = meanXX + current_base*current_base;
                        meanXY = meanXY + current_base*current_values;
                        counter = counter + 1;
                    }       
                }

                if(abs(current_max - current_min) < min_range){
                    output[y][x] = default_gradient;
                }

                else{
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