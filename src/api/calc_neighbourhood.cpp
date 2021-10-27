#include "gridpp.h"

#include <iostream>

using namespace gridpp;

vec2 gridpp::calc_neighbourhood(const vec2& array, const vec2& search_array,int halfwidth, 
    float search_criteria_min, float search_criteria_max , float search_target_min, float search_target_max, float search_delta){

    if(search_criteria_min > search_criteria_max){
        throw std::invalid_argument("Search_criteria_min must be smaller than search_criteria_max");
    }
    if(search_target_min > search_target_max){
        throw std::invalid_argument("Search_target_min must be smaller than search_target_max"); 
    }
    if(halfwidth < 0){
        throw std::invalid_argument("halfwidth must be positive");
    }

    vec2 output = gridpp::init_vec2(array.size(), array[0].size());

    int nY = array.size();
    int nX = array[0].size();

    
    for(int y = 0; y < array.size(); y++){ 
        for(int x = 0; x < array[y].size(); x++){
            /* Loop over each element in array */
            
            float nearest_target = gridpp::MV;
            int I_nearestSearchArray_Y = 0;
            int I_nearestSearchArray_X = 0;
            int counter = 0;
            float accum_temp = 0;

            if(!gridpp::is_valid(search_array[y][x])){
                /* if current search_array is invalid, set output equal to array (input) */
                output[y][x] = array[y][x];
                continue;
            }

            if(search_array[y][x] < search_criteria_min || search_array[y][x] >= search_criteria_max){
                /* Ignore values outside of range of scope (outside search_criteria min and max) */
                output[y][x] = array[y][x];
                continue;
            }


            for(int yy = std::max(0, y - halfwidth); yy <= std::min(nY - 1, y + halfwidth); yy++){
                for(int xx = std::max(0, x - halfwidth); xx <= std::min(nX - 1, x + halfwidth); xx++){
                    /* Loop over neighbourhood of y and x */
                    int I_SearchArray_Y = 0;
                    int I_SearchArray_X = 0;

                    if(!gridpp::is_valid(search_array[yy][xx]) || !gridpp::is_valid(array[yy][xx])){
                        continue;
                    }

                    else if(search_array[y][x] >= search_criteria_min && search_array[y][x] <= search_criteria_max){
                        /*condition search array inside of scope */

                        if(search_array[yy][xx] >= search_target_min && search_array[yy][xx] <= search_target_max){
                            /* Count all and add all values*/
                            counter++;
                            accum_temp = accum_temp + array[yy][xx];
                        }
                        else if (counter > 0)
                            /*  We have decided that we will only use values within th search target */
                            continue;
                        else if(std::abs(search_array[yy][xx] - search_array[y][x]) >= search_delta) {
                            /* Finding nearest value that's outside of the search target range */
                            if(!gridpp::is_valid(nearest_target)){
                                /* Set first value*/
                                nearest_target = search_array[yy][xx];
                                I_nearestSearchArray_Y = yy;
                                I_nearestSearchArray_X = xx;
                            }
                            else {
                                float curr_dist_to_target = std::min(std::abs(search_array[yy][xx] - search_target_min), std::abs(search_array[yy][xx] - search_target_max));
                                float best_dist_to_target = std::min(std::abs(nearest_target - search_target_min), std::abs(nearest_target - search_target_max));
                                if(curr_dist_to_target < best_dist_to_target) {
                                    // If next search array is closer to search target, assign new value*
                                    nearest_target = search_array[yy][xx];
                                    I_nearestSearchArray_Y = yy;
                                    I_nearestSearchArray_X = xx;
                                }
                            }
                        }
                    }
                }
            }

            if(counter > 0){
                /* find mean value of accumulated values */
                output[y][x] = accum_temp / counter;
            }

            else if(gridpp::is_valid(nearest_target)){
                /* If no values found inside target, and nearest target was used, assign the value from that location*/
                output[y][x] = array[I_nearestSearchArray_Y][I_nearestSearchArray_X];
            }

            else{
                /* If no methods use, use original value*/
                output[y][x] = array[y][x];
            }
        }
    }
    return output;
}
