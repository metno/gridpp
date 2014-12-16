#ifndef UTIL_H
#define UTIL_H
#include <string>

class Util {
   public:
      Util();
      static void error(std::string iMessage);
      static void warning(std::string iMessage);
};
#endif
