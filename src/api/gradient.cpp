#include "gridpp.h"

using namespace gridpp;

vec2 gridpp::full_gradient(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient, Downscaler downscaler) {

    // Sizes;
    int nY = ogrid.size()[0];
    int nX = ogrid.size()[1];


    if(ivalues.size() != igrid.size()[0] || ivalues[0].size() != igrid.size()[1])
        throw std::invalid_argument("Values is the wrong size");

    if(laf_gradient.size() > 0)
        if(laf_gradient.size() != igrid.size()[0] || laf_gradient[0].size() != igrid.size()[1])

            throw std::invalid_argument("Laf gradient is the wrong size");

    if(elev_gradient.size() > 0)
        if(elev_gradient.size() != igrid.size()[0] || elev_gradient[0].size() != igrid.size()[1])
            throw std::invalid_argument("Elevation gradient is the wrong size");

    //Outputs
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();

    vec2 output = gridpp::downscaling(igrid, ogrid, ivalues, downscaler);
    vec2 delev_gradient;
    vec2 delevs;
    if(elev_gradient.size() != 0) {
        delev_gradient = gridpp::downscaling(igrid, ogrid, elev_gradient, downscaler);
        delevs = gridpp::downscaling(igrid, ogrid, igrid.get_elevs(), downscaler);
    }
    vec2 dlaf_gradient;
    vec2 dlafs;
    if(laf_gradient.size() != 0) {
        dlaf_gradient = gridpp::downscaling(igrid, ogrid, laf_gradient, downscaler);
        dlafs = gridpp::downscaling(igrid, ogrid, igrid.get_lafs(), downscaler);
    }

    #pragma omp parallel for collapse(2)
    for(int y = 0; y < nY ; y++){
        for(int x = 0; x < nX; x++){
            //Calculate LAF and elevation difference between output and input
            float laf_correction = 0;
            if(dlaf_gradient.size() > 0) {
                float olaf = olafs[y][x];
                float ilaf = dlafs[y][x];
                if(gridpp::is_valid(olaf) && gridpp::is_valid(ilaf)) {
                    float laf_diff = olaf - ilaf;
                    laf_correction = dlaf_gradient[y][x]*laf_diff;
                }
            }

            float elev_correction = 0;
            if(delev_gradient.size() > 0) {
                float oelev = oelevs[y][x];
                float ielev = delevs[y][x];
                if(gridpp::is_valid(oelev) && gridpp::is_valid(ielev)) {
                    float elev_diff = oelev - ielev;
                    elev_correction = delev_gradient[y][x]*elev_diff;
                }
            }

            output[y][x] += laf_correction + elev_correction;
        }
    }
    return output;
}

vec3 gridpp::full_gradient(const Grid& igrid, const Grid& ogrid, const vec3& ivalues, const vec3& elev_gradient, const vec3& laf_gradient, Downscaler downscaler) {
    // Sizes;
    int nY = ogrid.size()[0];
    int nX = ogrid.size()[1];
    int nTime = ivalues.size();

    //Outputs
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();

    vec3 output = gridpp::downscaling(igrid, ogrid, ivalues, downscaler);
    vec2 delevs;
    vec3 delev_gradient;
    if(elev_gradient.size() != 0) {
        assert(gridpp::compatible_size(elev_gradient, ivalues));
        delev_gradient = gridpp::downscaling(igrid, ogrid, elev_gradient, downscaler);
        delevs = gridpp::downscaling(igrid, ogrid, igrid.get_elevs(), downscaler);
    }
    vec3 dlaf_gradient;
    vec2 dlafs;
    if(laf_gradient.size() != 0) {
        assert(gridpp::compatible_size(laf_gradient, ivalues));
        dlaf_gradient = gridpp::downscaling(igrid, ogrid, laf_gradient, downscaler);
        dlafs = gridpp::downscaling(igrid, ogrid, igrid.get_lafs(), downscaler);
    }

    #pragma omp parallel for collapse(2)
    for(int y = 0; y < nY ; y++){
        for(int x = 0; x < nX; x++){
            for(int t = 0; t < nTime; t++){
                float laf_correction = 0;
                if(dlaf_gradient.size() > 0) {
                    float olaf = olafs[y][x];
                    float ilaf = dlafs[y][x];
                    if(gridpp::is_valid(olaf) && gridpp::is_valid(ilaf)) {
                        float laf_diff = olaf - ilaf;
                        laf_correction = dlaf_gradient[t][y][x]*laf_diff;
                    }
                }

                float elev_correction = 0;
                if(delev_gradient.size() > 0) {
                    float oelev = oelevs[y][x];
                    float ielev = delevs[y][x];
                    if(gridpp::is_valid(oelev) && gridpp::is_valid(ielev)) {
                        float elev_diff = oelev - ielev;
                        elev_correction = delev_gradient[t][y][x]*elev_diff;
                    }
                }

                output[t][y][x] += laf_correction + elev_correction;
            }
        }
    }
    return output;
}

vec gridpp::full_gradient(const Grid& igrid, const Points& opoints, const vec2& ivalues, const vec2& elev_gradient, const vec2& laf_gradient, Downscaler downscaler) {
    vec olafs = opoints.get_lafs();
    vec oelevs = opoints.get_elevs();

    int nPoints = opoints.size();

    vec output = gridpp::downscaling(igrid, opoints, ivalues, downscaler);
    vec delev_gradient;
    vec delevs;
    if(elev_gradient.size() != 0) {
        assert(elev_gradient.size() == ivalues.size());
        assert(elev_gradient[0].size() == ivalues[0].size());
        delev_gradient = gridpp::downscaling(igrid, opoints, elev_gradient, downscaler);
        delevs = gridpp::downscaling(igrid, opoints, igrid.get_elevs(), downscaler);
    }
    vec dlaf_gradient;
    vec dlafs;
    if(laf_gradient.size() != 0) {
        assert(laf_gradient.size() == ivalues.size());
        assert(laf_gradient[0].size() == ivalues[0].size());
        dlaf_gradient = gridpp::downscaling(igrid, opoints, laf_gradient, downscaler);
        dlafs = gridpp::downscaling(igrid, opoints, igrid.get_lafs(), downscaler);
    }

    for(int i = 0; i < nPoints; i++){
        //Calculate LAF and elevation difference between output and input
        float laf_correction = 0;
        if(dlaf_gradient.size() > 0) {
            float olaf = olafs[i];
            float ilaf = dlafs[i];
            if(gridpp::is_valid(olaf) && gridpp::is_valid(ilaf)) {
                float laf_diff = olaf - ilaf;
                laf_correction = dlaf_gradient[i]*laf_diff;
            }
        }

        float elev_correction = 0;
        if(delev_gradient.size() > 0) {
            float oelev = oelevs[i];
            float ielev = delevs[i];
            if(gridpp::is_valid(oelev) && gridpp::is_valid(ielev)) {
                float elev_diff = oelev - ielev;
                elev_correction = delev_gradient[i]*elev_diff;
            }
        }

        output[i] += laf_correction + elev_correction;
    }
    return output;
}

vec2 gridpp::full_gradient(const Grid& igrid, const Points& opoints, const vec3& ivalues, const vec3& elev_gradient, const vec3& laf_gradient, Downscaler downscaler) {
    vec olafs = opoints.get_lafs();
    vec oelevs = opoints.get_elevs();

    int nPoints = opoints.size();
    int nTime = ivalues.size();

    vec2 output = gridpp::downscaling(igrid, opoints, ivalues, downscaler);
    vec2 delev_gradient;
    vec delevs;
    if(elev_gradient.size() != 0) {
        delev_gradient = gridpp::downscaling(igrid, opoints, elev_gradient, downscaler);
        delevs = gridpp::downscaling(igrid, opoints, igrid.get_elevs(), downscaler);
    }
    vec2 dlaf_gradient;
    vec dlafs;
    if(laf_gradient.size() != 0) {
        dlaf_gradient = gridpp::downscaling(igrid, opoints, laf_gradient, downscaler);
        dlafs = gridpp::downscaling(igrid, opoints, igrid.get_lafs(), downscaler);
    }

    for(int i = 0; i < nPoints; i++){
        for(int t = 0; t < nTime; t++){
            float laf_correction = 0;
            if(dlaf_gradient.size() > 0) {
                float olaf = olafs[i];
                float ilaf = dlafs[i];
                if(gridpp::is_valid(olaf) && gridpp::is_valid(ilaf)) {
                    float laf_diff = olaf - ilaf;
                    laf_correction = dlaf_gradient[t][i]*laf_diff;
                }
            }

            float elev_correction = 0;
            if(delev_gradient.size() > 0) {
                float oelev = oelevs[i];
                float ielev = delevs[i];
                if(gridpp::is_valid(oelev) && gridpp::is_valid(ielev)) {
                    float elev_diff = oelev - ielev;
                    elev_correction = delev_gradient[t][i]*elev_diff;
                }
            }

            output[t][i] += laf_correction + elev_correction;
        }
    }
    return output;
}

vec3 gridpp::full_gradient_debug(const Grid& igrid, const Grid& ogrid, const vec2& ivalues,  const vec2& elev_gradient, const vec2& laf_gradient, Downscaler downscaler) {
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
    vec2 ielevs = igrid.get_elevs();
    vec2 ilafs = igrid.get_lafs();

    vec2 i_lafs_gradient = laf_gradient;
    vec2 i_elevs_gradient = elev_gradient;

    //Outputs
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();

    vec3 output = gridpp::init_vec3(nY, nX, 9);
    vec2 dvalues = gridpp::downscaling(igrid, ogrid, ivalues, downscaler);
    vec2 delev_gradient;
    vec2 delevs;
    if(elev_gradient.size() != 0) {
        delev_gradient = gridpp::downscaling(igrid, ogrid, elev_gradient, downscaler);
        delevs = gridpp::downscaling(igrid, ogrid, igrid.get_elevs(), downscaler);
    }
    vec2 dlaf_gradient;
    vec2 dlafs;
    if(laf_gradient.size() != 0) {
        dlaf_gradient = gridpp::downscaling(igrid, ogrid, laf_gradient, downscaler);
        dlafs = gridpp::downscaling(igrid, ogrid, igrid.get_lafs(), downscaler);
    }

    for(int y = 0; y < nY ; y++){
        for(int x = 0; x < nX; x++){
            //Collect output and input LAF and elevations
            float olaf = olafs[y][x];
            float ilaf = dlafs[y][x];
            float oelev = oelevs[y][x];
            float ielev = delevs[y][x];

            //Calculate LAF and elevation difference between output and input
            float laf_diff = olaf - ilaf;
            float elev_diff = oelev - ielev;

            float laf_correction = 0;
            if(dlaf_gradient.size() > 0) {
                if(gridpp::is_valid(olaf) && gridpp::is_valid(ilaf)) {
                    float laf_diff = olaf - ilaf;
                    laf_correction = dlaf_gradient[y][x]*laf_diff;
                }
            }

            float elev_correction = 0;
            if(delev_gradient.size() > 0) {
                if(gridpp::is_valid(oelev) && gridpp::is_valid(ielev)) {
                    float elev_diff = oelev - ielev;
                    elev_correction = delev_gradient[y][x]*elev_diff;
                }
            }

            //Calculate temperature
            float temp = dvalues[y][x] + laf_correction + elev_correction;
            float temp_laf_c = dvalues[y][x] + laf_correction;
            float temp_elev_c = dvalues[y][x]+ elev_correction;

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
