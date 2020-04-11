#include "gridpp.h"
float gridpp::pressure(float ielev, float oelev, float ipressure, float itemperature) {
    float g0 = 9.80665;
    float M = 0.0289644;
    float R = 8.3144598;
    float value = gridpp::MV;
    if(gridpp::util::is_valid(ielev) && gridpp::util::is_valid(oelev) && gridpp::util::is_valid(ipressure) && gridpp::util::is_valid(itemperature))
        value = ipressure * exp(-g0 * M * (oelev-ielev) / (R * itemperature));
    return value;
}

