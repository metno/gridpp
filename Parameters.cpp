#include "Parameters.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include "Util.h"

Parameters::Parameters(std::vector<float> iValues) : mValues(iValues) {

}
Parameters::Parameters() {

}

float const& Parameters::operator[](unsigned int i) const {
   if(i >= mValues.size()) {
      std::stringstream ss;
      ss << "Attempted to access element " << i << " in parameter array of length " << mValues.size() << std::endl;
      Util::error(ss.str());
   }
      
   return mValues[i];
}
int Parameters::size() const {
   return mValues.size();
}
