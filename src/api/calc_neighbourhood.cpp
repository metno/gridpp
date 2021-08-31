#include "gridpp.h"

#include <iostream>

using namespace gridpp;

vec2 gridpp::calc_neighbourhood(const vec2& array, const vec2& search_array,int halfwidth, 
    float search_criteria_min, float search_criteria_max , float search_target_min, float search_target_max){

    std::cout << "starting";

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
            
            bool start = true;

            float current_max = 0;
            int I_maxSearchArray_Y = 0;
            int I_maxSearchArray_X = 0;

            if(search_array[y][x] <= search_criteria_min){
                output[y][x] = array[y][x];
                continue;
            }
            
            if(search_array[y][x] > search_target_max){
                output[y][x] = array[y][x];
                continue;
            }

            for(int yy = std::max(0, y - halfwidth); yy <= std::min(nY - 1, y + halfwidth); yy++){
                for(int xx = std::max(0, x - halfwidth); xx <= std::min(nX - 1, x + halfwidth); xx++){
                    
                    float current_array = search_array[yy][xx];
                    if(!gridpp::is_valid(current_array)){
                        continue;
                    } 

                    if(start){
                        start = false;
                        current_max = current_array;
                        I_maxSearchArray_Y = yy;
                        I_maxSearchArray_X = xx;
                    }

                    else if(current_array > current_max){
                        current_max = current_array;
                        I_maxSearchArray_Y = yy;
                        I_maxSearchArray_X = xx;
                    }
                    
                    //std::cout << "inner loop";
                    
                }
            }
            
            
            if(!gridpp::is_valid(current_max)){
                output[y][x] = 0;
            }

            else{
                output[y][x] = array[I_maxSearchArray_Y][I_maxSearchArray_X];
            }
        }
    }
    std::cout << " ends ";
    return output;
}
/*  array = input temperature
*        search_array = input_lafs
*       halfwidth = radius to search for
*      search_criteria = value in which search is instigated 
*     search_target = accpetance value of laf.
*/


 /* 
 array (vec2)
 laf ( vec2) (search array)
 halfwidth (int)
 max_laf (float) search criteria
 min_laf (float) search target
*/