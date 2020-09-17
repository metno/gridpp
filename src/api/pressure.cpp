#include "gridpp.h"

using namespace gridpp;

float gridpp::pressure(float ielev, float oelev, float ipressure, float itemperature) {
    float g0 = 9.80665;
    float M = 0.0289644;
    float R = 8.3144598;
    float value = gridpp::MV;
    if(gridpp::is_valid(ielev) && gridpp::is_valid(oelev) && gridpp::is_valid(ipressure) && gridpp::is_valid(itemperature))
        value = ipressure * exp(-g0 * M * (oelev-ielev) / (R * itemperature));
    return value;
}

