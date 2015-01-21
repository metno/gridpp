#include "Options.h"

Options::Options() {

}
void Options::addOption(std::string iKey, std::string iValue) {
   mKeys.push_back(iKey);
   mValues.push_back(iValue);
}

void Options::addOption(std::string iOptionString) {
   int nextEqual = iOptionString.find('=');
   if(nextEqual < 0)
      return;
   std::string key   = iOptionString.substr(0, nextEqual);
   std::string value = iOptionString.substr(nextEqual+1);

   addOption(key, value);
}
void Options::clear() {
   mKeys.clear();
   mValues.clear();
}
