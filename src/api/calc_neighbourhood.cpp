#include "gridpp.h"

using namespace gridpp;

vec2 gridpp::calc_neighbourhood(const vec2& array, const vec2& search_array,int halfwidth, int search_criteria , int search_target){


    if(halfwidth < 0){
        throw std::invalid_argument("halfwidth must be positive");
    }

    int nY = array.size();
    int nX = array[0].size();

    vec2 output = gridpp::init_vec2(array.size(), array[0].size());

    for(int y = 0; y < array.size(); y++){
        for(int x = 0; x < array[y].size(); x++){
            bool start = true;

            float current_max = 0;
            int I_maxSearchArray_Y = 0;
            int I_maxSearchArray_X = 0;

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