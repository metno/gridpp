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

vec gridpp::pressure(const vec& ielev, const vec& oelev, const vec& ipressure, const vec& itemperature) {
    int N = ielev.size();
    if(oelev.size() != N || ipressure.size() != N || itemperature.size() != N)
        throw std::invalid_argument("pressure: Input arguments must be of the same size");

    vec values(N);
    #pragma omp parallel for
    for(int i = 0; i < N; i++) {
        values[i] = gridpp::pressure(ielev[i], oelev[i], ipressure[i], itemperature[i]);
    }
    return values;
}

float gridpp::sea_level_pressure(float ps, float altitude, float temperature, float rh, float dewpoint) {
    // For altitude it is assumed geopotential meters is close to meters.
    // Following https://library.wmo.int/doc_num.php?explnum_id=10616

    if ( !gridpp::is_valid(altitude) ) {
        throw std::runtime_error("sea_level_pressure: altitude is NAN");
    } else if ( !gridpp::is_valid(temperature) ){
        throw std::runtime_error("sea_level_pressure: temperature is NAN");
    } else if ( ps < 0. || temperature < 0. || rh < 0. || rh > 1. || dewpoint < 0. ){
        throw std::runtime_error("sea_level_pressure: unphysical values in input");
    }
    
    

    float T = temperature - 273.15;
    float Ts = 273.15 + T;
    float g = gridpp::gravit;
    float R = gridpp::gas_constant_si;
    float e = 0.;
    float slp = 0.;
    float a = 0.0065; // [K/gpm]
    float Ch = 0.12; // [K/hPa]
    ps = ps*0.01;

    if ( gridpp::is_valid(rh) ) {
        float es = 6.11 * pow(10., ((7.5 * T) / (237.3 + T)) );
        e = rh * es;
        float A = 17.625;
        float B = 243.04;
        float C = 6.1094;
        // http://climate.envsci.rutgers.edu/pdf/LawrenceRHdewpointBAMS.pdf
        dewpoint = (B * log(e / C)) / (A - log(e / C));
    } else if ( gridpp::is_valid(dewpoint) ) {
        dewpoint = dewpoint - 273.15;
        e = 6.11 * pow(10, ((7.5 * dewpoint) / (237.3 + dewpoint)));
        // float es = 6.11 * pow(10, ((7.5 * T) / (237.3 + T)));
        // rh = e/es;
    } else if ( !gridpp::is_valid(rh) && !gridpp::is_valid(dewpoint) ) {
        dewpoint = T - 3.;
    }

    if (altitude >= 50.) {
        slp = ps * exp((g * altitude / R) / (Ts + 0.5 * a * altitude + e * Ch));
    } else if (altitude < 50.) {
        float Tv = (273.15 + T) / (1 - 0.379 * (6.11 * pow(10.,  (((7.5 * dewpoint) / (237.7 + dewpoint))) )  / ps ));
        float Ck = ps * altitude / (29.27 * Tv);
        slp = ps + Ck;
    }

    slp = slp*100.; // [hPa] -> [Pa]
    
    return slp;
}

vec gridpp::sea_level_pressure(const vec& ps, const vec& altitude, const vec&  temperature, const vec& rh, const vec& dewpoint) {
    int N = ps.size();
    if(altitude.size() != N || temperature.size() != N || rh.size() != N || dewpoint.size() != N)
        throw std::invalid_argument("slp: Input arguments must be of the same size");

    vec values(N);
    #pragma omp parallel for
    for(int i = 0; i < N; i++) {
        values[i] = gridpp::sea_level_pressure(ps[i], altitude[i], temperature[i], rh[i], dewpoint[i]);
    }
    return values;
}
