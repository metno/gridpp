#include "Util.h"
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

Util::Util() {

}
void Util::error(std::string iMessage) {
   std::cout << iMessage << std::endl;
   abort();
}
void Util::warning(std::string iMessage) {
   std::cout << iMessage << std::endl;
   abort();
}

double Util::clock() {
   timeval t;
   gettimeofday(&t, NULL);
   double sec = (t.tv_sec);
   double msec= (t.tv_usec);
   return sec + msec/1e6;
}
