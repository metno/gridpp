#ifndef FIELD_H
#define FIELD_H
#include <iostream>
#include <vector>

class Field {
   public:
      Field(std::vector<std::vector<float> > iValues);
   private:
      std::string mVariable;
      std::vector<std::vector<float> > mValues;
};
#endif

