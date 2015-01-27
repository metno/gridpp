#ifndef VARIABLE_H
#define VARIABLE_H
#include <iostream>
#include <vector>
class Variable {
   public:
      // Precip: Hourly amount ending at the specified time
      // PrecipAcc: Accumulated ending at the specified time
      enum Type {Precip = 0, PrecipAcc = 1, Cloud = 10, T = 20, U = 30, V = 40, W = 50};
      static std::string getTypeName(Type iType);
      static Type getType(std::string iName);
      static std::string description();
      static std::vector<Type> getAllVariables();
};
#endif
