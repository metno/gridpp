#include "Util.h"
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <assert.h>
bool Util::mShowError = false;
bool Util::mShowWarning = false;
bool Util::mShowStatus = false;
float Util::MV = -999;
float Util::pi  = 3.14159265;
double Util::radiusEarth = 6.378137e6;
Util::Util() {

}
void Util::error(std::string iMessage) {
   if(mShowError)
      std::cout << "ERROR:   " << iMessage << std::endl;
   abort();
}
void Util::warning(std::string iMessage) {
   if(mShowWarning)
      std::cout << "WARNING: " << iMessage << std::endl;
}
void Util::status(std::string iMessage) {
   if(mShowStatus)
      std::cout << "STATUS:  " << iMessage << std::endl;
}

double Util::clock() {
   timeval t;
   gettimeofday(&t, NULL);
   double sec = (t.tv_sec);
   double msec= (t.tv_usec);
   return sec + msec/1e6;
}

void Util::setShowError(bool flag) {
   mShowError = flag;
}

void Util::setShowWarning(bool flag) {
   mShowWarning = flag;
}

void Util::setShowStatus(bool flag) {
   mShowStatus = flag;
}

bool Util::isValid(float iValue) {
   return !std::isnan(iValue) && !std::isinf(iValue) && iValue != Util::MV;
}

float Util::getDistance(float lat1, float lon1, float lat2, float lon2) {
   if(!Util::isValid(lat1) || !Util::isValid(lat2) ||
      !Util::isValid(lon1) || !Util::isValid(lon2)) {
      return Util::MV;
   }
   assert(fabs(lat1) <= 90 && fabs(lat2) <= 90 && fabs(lon1) <= 360 && fabs(lon2) <= 360);

   if(lat1 == lat2 && lon1 == lon2)
      return 0;
   double lat1r = deg2rad(lat1);
   double lat2r = deg2rad(lat2);
   double lon1r = deg2rad(lon1);
   double lon2r = deg2rad(lon2);
   double ratio = cos(lat1r)*cos(lon1r)*cos(lat2r)*cos(lon2r) + cos(lat1r) * sin(lon1r) *cos(lat2r)*sin(lon2r) + sin(lat1r)*sin(lat2r);
   double dist = acos(ratio)*radiusEarth;
   return (float) dist;
}

float Util::deg2rad(float deg) {
   return (deg * Util::pi / 180);
}
float Util::rad2deg(float rad) {
   return (rad * 180 / Util::pi);
}
