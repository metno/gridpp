#ifndef VARIABLE_H
#define VARIABLE_H
#include <iostream>
class Variable {
   public:
      enum Type {Precip = 0, PrecipAcc = 1, Cloud = 10, T = 20};
      static std::string getTypeName(Type iType);
};
#endif
