#include "Options.h"

void Options::addOption(std::string iKey, std::string iValue) {
   if(iKey == "" || iValue == "")
      return;

   // Overwrite existing value if key already exists
   for(int i = 0; i < mKeys.size(); i++) {
      if(mKeys[i] == iKey) {
         mValues[i] = iValue;
         return;
      }
   }
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
