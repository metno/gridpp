#include "Variable.h"
#include <sstream>
#include "Util.h"

Variable::Variable() {
}
Variable::Variable(std::string iName, float iMin, float iMax, std::string iUnits, std::string iStandardName):
      mName(iName),
      mMin(iMin),
      mMax(iMax),
      mUnits(iUnits),
      mStandardName(iStandardName) {
}
Variable::Variable(std::string iName, std::string iUnits, std::string iStandardName):
      mName(iName),
      mMin(Util::MV),
      mMax(Util::MV),
      mUnits(iUnits),
      mStandardName(iStandardName) {
}

float Variable::min() const {
   return mMin;
}

float Variable::max() const {
   return mMax;
}

std::string Variable::name() const {
   return mName;
}

std::string Variable::units() const {
   return mUnits;
}

std::string Variable::standardName() const {
   return mStandardName;
}

bool Variable::operator<(const Variable &right) const {
   return mName < right.name();
}
