#ifndef VARIABLE_H
#define VARIABLE_H
#include <iostream>
#include <vector>

//! Represents a meteorological variable and its metadata
class Variable {
   public:

      Variable();
      Variable(std::string iName, float iMin, float iMax, std::string iUnits, std::string iStandardName);
      Variable(std::string iName, std::string iUnits, std::string iStandardName);

      // Name of variable
      std::string name() const;

      //! Get the minimum possible attainable value for this variable
      float min() const;

      //! Get the maximum possible attainable value for this variable
      float max() const;

      //! Returns the units of the variable
      std::string units() const;

      //! Returns the NetcdfCF standard name
      std::string standardName() const;

      bool operator<(const Variable &right) const;

   private:
      float mMin;
      float mMax;
      std::string mUnits;
      std::string mStandardName;
      std::string mName;
};
#endif
