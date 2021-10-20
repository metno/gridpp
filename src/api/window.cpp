#include "gridpp.h"
#include <iostream>

using namespace gridpp;


vec2 gridpp::window(const vec2& array, 
    int length, gridpp::Statistic statistic, bool before, 
    bool keep_missing, bool missing_edges){

    vec2 output = gridpp::init_vec2(array.size(), array[0].size(), 0);

    int nY = array.size();
    int nX = array[0].size();

    for(int y = 0; y < array.size(); y++){
        for(int x = 0; x < array[y].size(); x++){

            float accum = 0;
            float counter = 0;

            int start;
            int end;

            float max_value = gridpp::MV;
            float min_value = gridpp::MV;

            if(before == true){
                // BEFORE BOOLEAN: window end at the timestep
                start = x - length + 1;
                end = x;
            }
            else{
                // IF FALSE: centre the window on each timestep
                start = x - length / 2;
                end = x + length / 2;
            }

            for(int xx = start; xx <= end; xx++){
                if(xx < 0 || xx >= nX){
                    if(missing_edges == true){
                        // MISSING_EDGES BOOLEAN: set the window value to missing, 
                        // if the window goes outside the edges
                        accum = std::nanf("1");
                        break;
                    }
                    else{
                        continue;
                    }
                }
     
                if(std::isnan(array[y][xx])){
                    if(keep_missing == true){
                        // KEEP_MISSING BOOLEAN: set window value to missing if
                        // one or more values in window is missing
                        accum = std::nanf("1");
                        break;
                    }
                    else{
                        continue;
                    }
                }

                else{
                    accum = accum + array[y][xx];
                    counter++;

                    // MAX VALUE
                    if(!gridpp::is_valid(max_value)){
                        max_value = array[y][xx];
                    }
                    else if(xx > start){
                        if(array[y][xx] > array[y][xx-1]){
                            max_value = array[y][xx];
                        }           
                    }
                    // MIN VALUE
                    if(!gridpp::is_valid(min_value)){
                        min_value = array[y][xx];
                    }
                    else if(xx > start){
                        if(array[y][xx] < array[y][xx-1]){
                            min_value = array[y][xx];
                        }
                    }
                }
            }

            if(counter == 0){ 
                // If all values inside window are nans, output is nan.
                output[y][x] = std::nan("1");
                continue;
            }
            if(statistic == gridpp::Sum){
                // If statistic is Sum
                output[y][x] = accum;
            }  
            else if(statistic == gridpp::Mean){
                //If statistic is Mean
                output[y][x] = accum / counter;
            }
            else if(statistic == gridpp::Median){
                throw std::invalid_argument("Statistic Method not Implemented");
            }
            else if(statistic == gridpp::Max){
                output[y][x] = max_value;            
            }
            else if(statistic == gridpp::Min){
                output[y][x] = min_value;
            }
            else{
                throw std::invalid_argument("Statistic Method not implemented");
            }
        }
    }
    return output;
}