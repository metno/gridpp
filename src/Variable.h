#ifndef VARIABLE_H
#define VARIABLE_H
#include <iostream>
#include <vector>

// Represents a meteorological variable
class Variable {
   public:
      enum Type {
         Precip = 0,       // Hourly amount ending at the specified time
         PrecipAcc = 1,    // Accumulated ending at the specified time
         Cloud = 10,       // Cloud cover (between 0 and 1)
         T = 20,           // 2m temperature (K)
         U = 30,           // 10m U-wind (m/s)
         V = 40,           // 10m V-wind (m/s)
         W = 50,           // 10m windspeed (m/s)
         RH = 60,          // Reliative humidity (%)
         Phase = 70,       // Precip phase
         P = 80,           // Surface pressure (pa)
         Fake = 90         // Fake variable used for testing
      };
      //! Convert type to string
      static std::string getTypeName(Type iType);

      //! Get the minimum possible attainable value for this variable
      static float getMin(Type iType);

      //! Get the maximum possible attainable value for this variable
      static float getMax(Type iType);

      //! Convert string to type
      static Type getType(std::string iName);

      //! Description of all defined variables
      static std::string description();

      //! Vector of all defined variables
      static std::vector<Type> getAllVariables();
};
#endif
