#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <vector>

class Parameters {
   public:
      Parameters(std::vector<float> iValues);
      Parameters();
      int size() const;
      std::vector<float> getValues() const;
            float & operator[](unsigned int i);
      const float & operator[](unsigned int i) const;
   private:
      std::vector<float> mValues;
};
#endif

