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
            
            float current_max = gridpp::MV;

            int I_maxSearchArray_Y = 0;
            int I_maxSearchArray_X = 0;

            int counter = 0;
            float accum_temp = 0;

            if(search_array[y][x] < search_criteria_min || search_array[y][x] >= search_criteria_max){
                output[y][x] = array[y][x];
                continue;
            }

            for(int yy = std::max(0, y - halfwidth); yy <= std::min(nY - 1, y + halfwidth); yy++){
                for(int xx = std::max(0, x - halfwidth); xx <= std::min(nX - 1, x + halfwidth); xx++){
                    
                    float current_array = search_array[yy][xx];

                    int I_SearchArray_Y = 0;
                    int I_SearchArray_X = 0;

                    if(!gridpp::is_valid(current_array)){
                        continue;
                    } 

                    else if(!gridpp::is_valid(current_max)){
                        current_max = current_array;
                        I_maxSearchArray_Y = yy;
                        I_maxSearchArray_X = xx;
                    }

                    /*else if(current_array >= current_max){
                        current_max = current_array;
                        I_maxSearchArray_Y = yy;
                        I_maxSearchArray_X = xx;
                    }*/
                    else if(search_array[y][x] >= search_criteria_min && search_array[y][x] <= search_criteria_max){
                        if(current_array >= search_target_min && current_array <= search_target_max){
                            //accum = accum + current_array;
                            //I_SearchArray_Y = yy;
                            //I_SearchArray_X = xx;
                        
                            counter++; //counter = counter + 1
                            accum_temp = accum_temp + array[yy][xx];
                        }

                        if(current_array >= current_max){
                            //should this be else if??? consequence: current_max is maximum value outside the search_target
                            current_max = current_array;
                            I_maxSearchArray_Y = yy;
                            I_maxSearchArray_X = xx;
                        }
                    }

                    else{
                        continue;
                    }
                }
            }         
            
            if(!gridpp::is_valid(current_max)){
                output[y][x] = array[y][x];
            }

            /*else if(current_max < search_target_min){
                output[y][x] = array[y][x];
            }*/ 

            /*else{
                output[y][x] = array[I_maxSearchArray_Y][I_maxSearchArray_X];
            }*/

            else if(counter > 0){
                output[y][x] = accum_temp / counter;  
            }

            //else if(current_max < search_target_min){
            //    output[y][x] = array[I_maxSearchArray_Y][I_maxSearchArray_X];
            //} 

            else{
                output[y][x] = array[y][x];
            }
        }
    }
    std::cout << " ends ";
    return output;
}
