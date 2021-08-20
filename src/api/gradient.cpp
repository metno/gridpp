#include "gridpp.h"

using namespace gridpp;



vec2 gridpp::gradient(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient){
    //Inputs
    //vec2 ivalues; 
    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 ielevs = igrid.get_elevs();
    vec2 ilafs = igrid.get_lafs();

    vec2 i_lafs_gradient = laf_gradient;
    vec2 i_elevs_gradient = elev_gradient;

    //Outputs
    //vec2 o_values = ogrid::nearest(igrid, ogrid, ivalues); Dont need this?
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();



    // Sizes;
    //int nY = ogrid.size();
    //int nX = ogrid[0].size(); 

    int nY = ogrid.size()[0]; 
    int nX = ogrid.size()[1];

    vec2 output = gridpp::init_vec2(nY, nX);

    for(int y = 0; y < nY ; y++){
        for(int x = 0; x < nX; x++){
            //collect index for nearest neighbour
            ivec indices = igrid.get_nearest_neighbour(olats[y][x], olons[y][x]);

            //Collect output and input LAF and elevations
            float olaf = olafs[y][x];
            float ilaf = ilafs[indices[0]][indices[1]];
            float oelev = oelevs[y][x];
            float ielev = ielevs[indices[0]][indices[1]];

            //Calculate LAF and elevation difference between output and input
            float laf_diff = olaf - ilaf; 
            float elev_diff = oelev - ielev;

            float laf_correction = i_lafs_gradient[indices[0]][indices[1]]*laf_diff;
            float elev_correction = i_elevs_gradient[indices[0]][indices[1]]*elev_diff;

            //Calculate temperature
            float temp = ivalues[indices[0]][indices[1]] + laf_correction + elev_correction;

            //Assign value to output
            output[y][x] = temp;
        }
    }
    return output;
}


vec3 gridpp::correction(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient){
    //Inputs
    //vec2 ivalues; 

    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 ielevs = igrid.get_elevs();
    vec2 ilafs = igrid.get_lafs();

    vec2 i_lafs_gradient = laf_gradient;
    vec2 i_elevs_gradient = elev_gradient;

    //Outputs
    //vec2 o_values = ogrid::nearest(igrid, ogrid, ivalues); Dont need this?
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();

    // Sizes;
    int nY = ogrid.size()[0];
    int nX = ogrid.size()[1];

    vec3 output = gridpp::init_vec3(nY, nX, 9);

    for(int y = 0; y < nY ; y++){
        for(int x = 0; x < nX; x++){
            //collect index for nearest neighbour
            ivec indices = igrid.get_nearest_neighbour(olats[y][x], olons[y][x]);

            //Collect output and input LAF and elevations
            float olaf = olafs[y][x];
            float ilaf = ilafs[indices[0]][indices[1]];
            //ilaf = ilaf/100.0;
            float oelev = oelevs[y][x];
            float ielev = ielevs[indices[0]][indices[1]];

            //Calculate LAF and elevation difference between output and input
            float laf_diff = olaf - ilaf;
            float elev_diff = oelev - ielev;
            
            float laf_correction = i_lafs_gradient[indices[0]][indices[1]]*laf_diff;
            float elev_correction = i_elevs_gradient[indices[0]][indices[1]]*elev_diff;

            //Calculate temperature
            float temp = ivalues[indices[0]][indices[1]] + laf_correction + elev_correction;
            float temp_laf_c = ivalues[indices[0]][indices[1]] + laf_correction;
            float temp_elev_c = ivalues[indices[0]][indices[1]]+ elev_correction;


            //Assign value to output
            output[y][x][0] = temp;
            output[y][x][1] = laf_diff;
            output[y][x][2] = elev_diff;
            output[y][x][3] = laf_correction;
            output[y][x][4] = elev_correction;
            output[y][x][5] = olaf;
            output[y][x][6] = ilaf;
            output[y][x][7] = temp_laf_c;
            output[y][x][8] = temp_elev_c;                        
        }
    }
    return output;
}