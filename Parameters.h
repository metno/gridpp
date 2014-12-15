#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <vector>

class Parameters {
   public:
      Parameters(std::vector<float> iValues);
      std::vector<float> getValues() const;
   private:
      std::vector<float> mValues;
};
#endif

