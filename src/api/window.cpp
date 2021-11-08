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
        
        // Mean and Sum Statistic 
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

                else{
                    start = std::max(0, x - length / 2);
                    end = std::min(nX - 1, x + length / 2);
                }

                // 
                if(start - 1 >= 0){
                    if(counts[end] - counts[start-1] == 0){
                        output[y][x] = gridpp::MV;
                    }

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


        /*else if(statistic == gridpp:: Max || statistic == gridpp::Min){
            throw std::invalid_argument("Statistic currently not supported");
        }*/



        // METHOD: use built in calc_Statistic 
        /*
        else if(statistic == gridpp::Max || statistic == gridpp::Min){
            for(int x = 0; x < array[y].size(); x++){
                int start;
                int end;

                float minimum = gridpp::MV;
                float maximum = gridpp::MV;

                if(before == true){
                    start = std::max(0, x - length + 1);
                    end = x;
                }
                else{
                    start = std::max(0, x - length / 2);
                    end = std::min(nX - 1, x + length / 2);
                }

                std::vector <float> stat_array (end - start, 0);
                for(int i = 0; i <= end - start; i++){
                    if(0 > x - length / 2){
                        stat_array[i] = array[y][x + i];
                    }
                    else{
                        stat_array[i] = array[y][x + i - 1];

                    }
                }                

                std::cout << " --------- " << "\n";
                int array_size = (sizeof(stat_array));
                std::cout << array_size  << "\n";

                for(int i = 0; i <= end - start; i++){
                    if( 0 > x - length / 2){
                        std::cout << y << " " << x + i << " " << stat_array[i] << "\n";
                    }

                    else{
                        std::cout << y << " " << x + i - 1 << " " << stat_array[i] << "\n";
                    }
                }

                std::cout << " ===> ";
                std::cout << "\n";

                //if(keep_missing == true){
                //    std::vector<int>::iterator it;
                //    it = std::find(stat_array.begin(), stat_array.end(), gridpp::MV);
                //    if (it != stat_array.end()) {
                //        output[y][x] = gridpp::MV;
                //    }
                //}

                if(statistic == gridpp::Min){
                    minimum = gridpp::calc_statistic(stat_array, statistic);
                    output[y][x] = minimum;
                    std::cout << minimum << " " << output[y][x]; 
                }

                std::cout << "\n";
                std::cout << " :::::: ";
                std::cout << "\n";

                if(statistic == gridpp::Max){
                    maximum = gridpp::calc_statistic(stat_array, statistic);
                    output[y][x] = maximum;
                    std::cout << maximum;
                }


                if(missing_edges == true){
                    if(0 > x - length / 2 || nX - 1 < x + length / 2){
                        output[y][x] = gridpp::MV;
                    }
                }
            }
        }
        */



        
        // Max and Min statistic
        // Use Method: maunally calculate max and min 
        else if(statistic == gridpp::Max || statistic == gridpp::Min){
            for(int x = 0; x < array[y].size(); x++){

                float accum = 0;
                float counter = 0;

                int start;
                int end;

                float max_value = gridpp::MV;
                float min_value = gridpp::MV;

                float minimum = gridpp::MV;
                float maximum = gridpp::MV;

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

                // gridpp::calc_statistic(array[y][start:end], gridpp::Max)
                // gridpp::calc_statistic(array[y][start:end], gridpp::Min)

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
                else if(statistic == gridpp::Max){
                    //output[y][x] = max_value;     
                    output[y][x] = maximum;       
                }
                else if(statistic == gridpp::Min){
                    //output[y][x] = min_value;
                    output[y][x] = minimum; 
                }
            }            
        }
        
        else{
            throw std::invalid_argument("Statistic currently not supported");
        }
    }
    return output;
}