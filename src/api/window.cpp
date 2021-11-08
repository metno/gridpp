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
        else {
            for(int x = 0; x < array[y].size(); x++){
                int start;
                int end;
                bool outside = false;

                if(before == true){
                    start = x - length + 1;
                    end = x;
                }
                else {
                    start = x - length / 2;
                    end = x + length / 2;
                }
                if(start < 0) {
                    start = 0;
                    outside = true;
                }
                if(end > nX - 1) {
                    end = nX - 1;
                    outside = true;
                }
                int array_size = end - start + 1;

                vec stat_array(array_size, 0);
                int count_missing = 0;
                for(int i = start; i <= end; i++){
                    float curr = array[y][i];
                    count_missing += !gridpp::is_valid(curr);
                    stat_array[i - start] = curr;
                }

                //if(keep_missing == true){
                //    std::vector<int>::iterator it;
                //    it = std::find(stat_array.begin(), stat_array.end(), gridpp::MV);
                //    if (it != stat_array.end()) {
                //        output[y][x] = gridpp::MV;
                //    }
                //}

                if(keep_missing && count_missing > 0)
                    output[y][x] = gridpp::MV;
                else if(missing_edges && outside)
                    output[y][x] = gridpp::MV;
                else {
                    float value = gridpp::calc_statistic(stat_array, statistic);
                    output[y][x] = value;
                }
            }
        }
    }
    return output;
}
