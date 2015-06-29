#include "Parameters.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include "Util.h"

Parameters::Parameters(std::vector<float> iValues) : mValues(iValues) {

}
Parameters::Parameters(float iValue) {
   mValues = std::vector<float>(1, iValue);
}
Parameters::Parameters() {

}

std::vector<float> Parameters::getValues() const {
   return mValues;
}
float & Parameters::operator[](unsigned int i) {
   if(i >= mValues.size() || i < 0) {
      std::stringstream ss;
      ss << "Attempted to access element " << i << " in parameter array of length " << mValues.size() << std::endl;
      Util::error(ss.str());
   }
      
   return mValues[i];
}
const float & Parameters::operator[](unsigned int i) const {
   if(i >= mValues.size() || i < 0) {
      std::stringstream ss;
      ss << "Attempted to access element " << i << " in parameter array of length " << mValues.size() << std::endl;
      Util::error(ss.str());
   }
      
   return mValues[i];
}
int Parameters::size() const {
   return mValues.size();
}

bool Parameters::isValid() const {
   for(int i = 0; i < mValues.size(); i++) {
      if(!Util::isValid(mValues[i]))
         return false;
   }
   return true;
}
