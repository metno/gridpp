#include "gridpp.h"

using namespace gridpp;

vec2 gridpp::full_gradient(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient){

    // Sizes;
    int nY = ogrid.size()[0];
    int nX = ogrid.size()[1];

    //if(ivalues.size() != nY || ivalues[0].size() != nX)
    //    throw std::invalid_argument("Values is the wrong size");

    /*if(laf_gradient.size() > 0)
        if(laf_gradient.size() != nY || laf_gradient[0].size() != nX)
            throw std::invalid_argument("Laf gradient is the wrong size");

    if(elev_gradient.size() > 0)
        if(elev_gradient.size() != nY || elev_gradient[0].size() != nX)
            throw std::invalid_argument("Elevation gradient is the wrong size");
    */
    //Inputs
    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 ielevs = igrid.get_elevs();
    vec2 ilafs = igrid.get_lafs();

    //Outputs
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();

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

            float laf_correction = 0;
            float elev_correction = 0;
            if(laf_gradient.size() > 0)
                laf_correction = laf_gradient[indices[0]][indices[1]]*laf_diff;
            if(elev_gradient.size() > 0)
                elev_correction = elev_gradient[indices[0]][indices[1]]*elev_diff;

            //Calculate temperature
            float temp = ivalues[indices[0]][indices[1]] + laf_correction + elev_correction;

            //Assign value to output
            output[y][x] = temp;
        }
    }
    return output;
}

vec3 gridpp::full_gradient(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec3& elev_gradient, const vec3& laf_gradient){
    // Sizes;
    int nY = ogrid.size()[0];
    int nX = ogrid.size()[1];
    int nTime = ivalues.size();

    //Inputs
    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 ielevs = igrid.get_elevs();
    vec2 ilafs = igrid.get_lafs();

    //Outputs
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();

    vec3 output = gridpp::init_vec3(nTime, nY, nX); 

    for(int y = 0; y < nY ; y++){
        for(int x = 0; x < nX; x++){
            ivec indices = igrid.get_nearest_neighbour(olats[y][x], olons[y][x]);
            int I = indices[0];
            int J = indices[1];

            //Collect output and input LAF and elevations
            float olaf = olafs[y][x];
            float ilaf = ilafs[I][J];
            float oelev = oelevs[y][x];
            float ielev = ielevs[I][J];

            //Calculate LAF and elevation difference between output and input
            float laf_diff = olaf - ilaf;
            float elev_diff = oelev - ielev;


            for(int t = 0; t < nTime; t++){

                float laf_correction = 0;
                float elev_correction = 0;

                if(laf_gradient.size() > 0)
                    laf_correction = laf_gradient[t][I][J]*laf_diff;
                if(elev_gradient.size() > 0)
                    elev_correction = elev_gradient[t][I][J]*elev_diff;

                float temp = ivalues[t][I][J] + laf_correction + elev_correction;
                
                output[t][y][x] = temp;
            }
        }
    }
    return output;
}

vec gridpp::full_gradient(const Grid& igrid, const Points& opoints, const vec2& ivalues, const vec2& elev_gradient, const vec2& laf_gradient){
    vec olats = opoints.get_lats();
    vec olons = opoints.get_lons();
    vec olafs = opoints.get_lafs();
    vec oelevs = opoints.get_elevs();

    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 ilafs = igrid.get_lafs();
    vec2 ielevs = igrid.get_elevs();

    int nPoints = olats.size();

    vec output(nPoints);

    for(int i = 0; i < nPoints; i++){
        ivec indices = igrid.get_nearest_neighbour(olats[i], olons[i]);
        int I = indices[0];
        int J = indices[1];

        float olaf = olafs[i];
        float ilaf = ilafs[I][J];
        float oelev = oelevs[i];
        float ielev = ielevs[I][J];

        //Calculate LAF and elevation difference between output and input
        float laf_diff = olaf - ilaf;
        float elev_diff = oelev - ielev;

        float laf_correction = 0;
        float elev_correction = 0;

        if(laf_gradient.size() > 0)
            laf_correction = laf_gradient[I][J]*laf_diff;
        if(elev_gradient.size() > 0)    
            elev_correction = elev_gradient[I][J]*elev_diff;

        float temp = ivalues[I][J] + laf_correction + elev_correction;

        output[i] = temp;    
    }
    return output;
}

vec2 gridpp::full_gradient(const Grid& igrid, const Points& opoints, const vec3& ivalues, const vec3& elev_gradient, const vec3& laf_gradient){
    vec olats = opoints.get_lats();
    vec olons = opoints.get_lons();
    vec olafs = opoints.get_lafs();
    vec oelevs = opoints.get_elevs();

    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 ilafs = igrid.get_lafs();
    vec2 ielevs = igrid.get_elevs();

    int nPoints = olats.size();
    int nTime = ivalues.size();

    vec2 output = gridpp::init_vec2(nTime, nPoints);

    for(int i = 0; i < nPoints; i++){
        ivec indices = igrid.get_nearest_neighbour(olats[i], olons[i]);
        int I = indices[0];
        int J = indices[1];

        float olaf = olafs[i];
        float ilaf = ilafs[I][J];
        float oelev = oelevs[i];
        float ielev = ielevs[I][J];

        //Calculate LAF and elevation difference between output and input
        float laf_diff = olaf - ilaf;
        float elev_diff = oelev - ielev;

        for(int t; t <nTime; t++){
            float laf_correction = 0;
            float elev_correction = 0;

            if(laf_gradient.size() > 0)
                laf_correction = laf_gradient[t][I][J] * laf_diff;
            if(elev_gradient.size() > 0)    
                elev_correction = elev_gradient[t][I][J] * elev_diff;


            float temp = ivalues[t][I][J] + laf_correction + elev_correction;

            output[t][i] = temp;
        }
    }
    return output;
}

vec3 gridpp::full_gradient_debug(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient){
    // Sizes;
    int nY = ogrid.size()[0];
    int nX = ogrid.size()[1];

    /*
    if(ivalues.size() != nY || ivalues[0].size() != nX)
        throw std::invalid_argument("Values is the wrong size");

    if(laf_gradient.size() > 0)
        if(laf_gradient.size() != nY || laf_gradient[0].size() != nX)
            throw std::invalid_argument("Laf gradient is the wrong size");

    if(elev_gradient.size() > 0)
        if(elev_gradient.size() != nY || elev_gradient[0].size() != nX)
            throw std::invalid_argument("Elevation gradient is the wrong size");
    */

    //Inputs
    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 ielevs = igrid.get_elevs();
    vec2 ilafs = igrid.get_lafs();

    vec2 i_lafs_gradient = laf_gradient;
    vec2 i_elevs_gradient = elev_gradient;

    //Outputs
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();

    vec3 output = gridpp::init_vec3(nY, nX, 9);

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

            float laf_correction = 0;
            float elev_correction = 0;
            if(laf_gradient.size() > 0)
                laf_correction = laf_gradient[indices[0]][indices[1]]*laf_diff;
            if(elev_gradient.size() > 0)
                elev_correction = elev_gradient[indices[0]][indices[1]]*elev_diff;

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
