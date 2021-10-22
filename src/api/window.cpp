#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec2 gridpp::window(const vec2& array,
    int length, gridpp::Statistic statistic, bool before,
    bool keep_missing, bool missing_edges){

    vec2 output = gridpp::init_vec2(array.size(), array[0].size(), 0);

    int nY = array.size();
    int nX = array[0].size();

    if(length % 2 == 0 && before == false){
        throw std::invalid_argument("Length variable must be an odd number");
    }

    #pragma omp parallel for
    for(int y = 0; y < array.size(); y++){
        if (statistic == gridpp::Mean || statistic == gridpp::Sum) {
            vec values = vec(array[y].size(), 0);
            vec counts = vec(array[y].size(), 0);

            for(int x = 0; x < array[y].size(); x++){
                if(x == 0){
                    if(gridpp::is_valid(array[y][x])){
                        values[x] = array[y][x];
                        counts[x] = 1;
                    }
                    else{
                        values[x] = 0;
                        counts[x] = 0;
                    }
                }

                else{
                    if(gridpp::is_valid(array[y][x])){
                        values[x] = array[y][x] + values[x-1];
                        counts[x] = counts[x-1] + 1;
                    }
                    else{
                        values[x] = values[x-1];
                        counts[x] = counts[x-1];
                    }
                }
            }

            for(int x = 0; x < array[y].size(); x++){
                int start;
                int end;

                // Compute Start and End points
                if(before == true){
                    start = std::max(0, x - length + 1);
                    end = x;
                }
                // Implement later (even numbers)
                //else if(length % 2 == 0){ // && before == false
                //    start = std::max(0, x - length / 2 - 1);
                //    end = std::min(nx - 1, x + length / 2 - 1)
                //}
                else{
                    start = std::max(0, x - length / 2);
                    end = std::min(nX - 1, x + length / 2);
                }

                // 
                if(start - 1 >= 0){
                    if(counts[end] - counts[start-1] == 0){
                        output[y][x] = gridpp::MV;
                    }
                    // Implement later (even numbers)
                    //else if(length % 2 == 0 && before == false){
                    //    output[y][x] = (values[end] + values[start]) / 2 + (values[end - 1] - values[start]);
                    //}

                    else{
                        output[y][x] = values[end] - values[start - 1];
                    }
                }
                else{
                    if(counts[end] == 0){
                        output[y][x] = gridpp::MV;
                    }
                    else{
                        output[y][x] = values[end];
                    }
                }

                if(statistic == gridpp::Mean){
                    if(counts[end] == 0){
                        output[y][x] = gridpp::MV;
                    }
                    else{
                        output[y][x] = output[y][x] / (counts[end] - counts[start-1]);
                    }
                }

                if(keep_missing == true){
                    if(counts[end] - counts[start - 1] < (end - (start - 1 ))){
                        output[y][x] = gridpp::MV;
                    }
                }

                if(missing_edges == true){
                    if(before == false){
                        if(x < length / 2 || x + length / 2  + 1 > array[y].size()){
                            output[y][x] = gridpp::MV;
                        }
                    }
                    else{
                        if(x < length - 1){
                            output[y][x] = gridpp::MV;
                        }
                    }
                }
            }
        }
        else{
            throw std::invalid_argument("Statistic currently not supported");
        }
    }
    return output;
}


namespace {
    vec2 window_old(const vec2& array,
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
                            accum = gridpp::MV;
                            break;
                        }
                        else{
                            continue;
                        }
                    }
         
                    if(!gridpp::is_valid(array[y][xx])){
                        if(keep_missing == true){
                            // KEEP_MISSING BOOLEAN: set window value to missing if
                            // one or more values in window is missing
                            accum = gridpp::MV;
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
                    output[y][x] = gridpp::MV;
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
}
