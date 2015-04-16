#include "Util.h"
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <assert.h>
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/date_time/date_duration.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <execinfo.h>
#include <signal.h>
#include <fstream>
#include <istream>
#include <iomanip>
#ifdef DEBUG
extern "C" void __gcov_flush();
#endif
bool Util::mShowError = true;
bool Util::mShowWarning = false;
bool Util::mShowStatus = false;
float Util::MV = -999;
float Util::pi  = 3.14159265;
double Util::radiusEarth = 6.378137e6;

void Util::error(std::string iMessage) {
   if(mShowError) {
      std::cout << "Error: " << iMessage << std::endl;
      void *array[10];
      size_t size = backtrace(array, 10);
      std::cout << "Stack trace:" << std::endl;
      backtrace_symbols_fd(array, size, 2);
   }
#ifdef DEBUG
   __gcov_flush();
#endif
   abort();
}
void Util::warning(std::string iMessage) {
   if(mShowWarning)
      std::cout << "Warning: " << iMessage << std::endl;
}
void Util::status(std::string iMessage) {
   if(mShowStatus)
      std::cout << "Status: " << iMessage << std::endl;
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

bool Util::exists(const std::string& iFilename) {
   std::ifstream infile(iFilename.c_str());
   return infile.good();
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
   double ratio = cos(lat1r)*cos(lon1r)*cos(lat2r)*cos(lon2r)
                + cos(lat1r)*sin(lon1r)*cos(lat2r)*sin(lon2r)
                + sin(lat1r)*sin(lat2r);
   double dist = acos(ratio)*radiusEarth;
   return (float) dist;
}

float Util::deg2rad(float deg) {
   return (deg * Util::pi / 180);
}
float Util::rad2deg(float rad) {
   return (rad * 180 / Util::pi);
}
int Util::getDate(time_t iUnixTime) {
   boost::gregorian::date epoch(1970,1,1);
   boost::gregorian::date_duration diff(iUnixTime/86400);
   boost::gregorian::date newDate = epoch + diff;

   return newDate.year() * 10000 + newDate.month() * 100 + newDate.day();
}
int Util::getCurrentDate() {
    boost::gregorian::date today = boost::gregorian::day_clock::local_day();
    return today.year()*10000 + today.month()*100 + today.day();
}

std::string Util::getCurrentDateString() {
    boost::gregorian::date today = boost::gregorian::day_clock::local_day();
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(4) << today.year() << "-"
       << std::setfill('0') << std::setw(2) << today.month() << "-"
       << std::setfill('0') << std::setw(2) << today.day();
    return ss.str();
}

int Util::calcDate(int iDate, float iAddHours) {
   int year  = floor(iDate / 10000);
   int month = floor(iDate % 10000)/100;
   int day   = floor(iDate % 100);

   int offDay = floor((iAddHours)/24);

   boost::gregorian::date currDate(year, month, day);
   boost::gregorian::date_duration diff(offDay);
   boost::gregorian::date newDate = currDate + diff; 
   int newYear  = newDate.year();
   int newMonth = newDate.month();
   int newDay   = newDate.day();
   int returnDate = newYear * 10000 + newMonth * 100 + newDay;
   return returnDate;
}

std::vector<std::string> Util::split(std::string iString) {
   std::vector<std::string> strings;
   std::stringstream ss(iString);
   std::string currString;
   while(ss >> currString) {
      strings.push_back(currString);
   }
   return strings;
}
float Util::logit(float p) {
   if(!Util::isValid(p) || p <= 0 || p >= 1)
      return Util::MV;
   return log(p/(1-p));
}
float Util::invLogit(float x) {
   if(!Util::isValid(x))
      return Util::MV;
   return exp(x)/(exp(x)+1);
}

bool Util::copy(std::string iFrom, std::string iTo) {
   std::ifstream source(iFrom.c_str(), std::ios::binary);
   std::ofstream dest(iTo.c_str(), std::ios::binary);
   if(!source.good() || !dest.good())
      return false;

   dest << source.rdbuf();

   source.close();
   dest.close();
   return true;
}
