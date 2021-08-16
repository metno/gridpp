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
    vec2 olafs = ogrid.get_lats();

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
            float laf_diff = (olaf - ilaf) / 100.0;
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


vec2 gridpp::correction(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient){
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
    vec2 olafs = ogrid.get_lats();

    // Sizes;
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
            float laf_diff = (olaf - ilaf) / 100.0;
            float elev_diff = oelev - ielev;

            float laf_correction = i_lafs_gradient[indices[0]][indices[1]]*laf_diff;
            float elev_correction = i_elevs_gradient[indices[0]][indices[1]]*elev_diff;

            //Calculate temperature
            float temp = ivalues[indices[0]][indices[1]] + laf_correction + elev_correction;

            //Assign value to output
            output[y][x] = laf_correction;
        }
    }
    return output;
}