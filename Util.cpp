#include "Util.h"
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
bool Util::mShowError = false;
bool Util::mShowWarning = false;
bool Util::mShowStatus = false;
float Util::MV = -999;
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
   return iValue != Util::MV;
}
