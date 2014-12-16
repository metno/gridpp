#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <vector>

class Parameters {
   public:
      Parameters(std::vector<float> iValues);
      Parameters();
      std::vector<float> getValues() const;
      float const& operator[](unsigned int i) const;
   private:
      std::vector<float> mValues;
};
#endif

