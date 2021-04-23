#include "gridpp.h"

using namespace gridpp;

vec2 gridpp::calc_gradient(const vec2& base, const vec2& values, int halfwidth , int min_num, float min_range, float default_gradient) {
    // calculate maximum element and minimum element within halfwidth of each [i][j]
    int nY = base.size();
    int nX = base[0].size();

    vec2 output = gridpp::init_vec2(base.size(), base[0].size());

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

    		float diffBase = current_max - current_min;
    		float diffValues = values[I_maxBase_Y][I_maxBase_X] - values[I_minBase_Y][I_minBase_X];
    		output[y][x] = diffValues / diffBase;
    		}
    	}
    }
    return output;
}