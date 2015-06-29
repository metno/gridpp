#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <vector>

//! Represents a vector of parameter values
class Parameters {
   public:
      Parameters(std::vector<float> iValues);
      Parameters(float iValue);

      //! Initialize an empty set of size 0
      Parameters();

      //! Returns the number of parameters
      int size() const;
      std::vector<float> getValues() const;

      //! Are all parameters in set valid numbers?
      bool isValid() const;

      //! Access the i'th parameter
            float & operator[](unsigned int i);
      const float & operator[](unsigned int i) const;
   private:
      std::vector<float> mValues;
};
#endif

