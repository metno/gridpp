#include "gridpp.h"
#include <math.h>
float gridpp::qnh(float pressure, float altitude) {
   if(pressure == 0)
         return 0;
   else if(gridpp::util::is_valid(altitude) && gridpp::util::is_valid(pressure)) {
      float g  = 9.80665;   // m/s2
      float T0 = 288.15;    // K
      float L  = 0.0065;    // K/m
      // Method 1:
      // float dElev = 0 - altitude;
      // float M  = 0.0289644; // kg/mol
      // float R  = 8.31447;   // J/(mol•K)
      // float cp = 1007;      // J/(kg•K)
      // float constant = -g*M/R/T0;
      // float qnh = pressure * pow(1 - L*dElev/T0, g*M/R/L);

      // Method 2: http://www.hochwarth.com/misc/AviationCalculator.html
      float CRGas = 287.053; // [m^2/(s^2*K)] = [J/(kg*K)]
      float p0    = 101325;  // pa
      float qnh   = p0*pow(pow((pressure/p0), (CRGas*L)/g) + (altitude*L)/T0, g/(CRGas*L));
      return qnh;
   }
   else {
      return MV;
   }
}
vec gridpp::qnh(const vec& pressure, const vec& altitude) {
    gridpp::util::check_equal_size(pressure, altitude);

    vec ret(pressure.size());
    for(int y = 0; y < pressure.size(); y++) {
        ret[y] = gridpp::qnh(pressure[y], altitude[y]);
    }
    return ret;
}
