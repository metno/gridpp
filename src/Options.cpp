#include "Options.h"
#include "Util.h"

void Options::addOption(std::string iKey, std::string iValue) {
   if(iKey == "" || iValue == "")
      return;

   // Overwrite existing value if key already exists
   for(int i = 0; i < mPairs.size(); i++) {
      std::string key = mPairs[i].first;
      if(key == iKey) {
         mPairs[i].second = iValue;
         return;
      }
   }
   mPairs.push_back(KeyValue(iKey, iValue));
}

void Options::addOptions(std::string iOptionString) {
   std::vector<std::string> options = Util::split(iOptionString);
   for(int i = 0; i < options.size(); i++) {
      addOption(options[i]);
   }
}
void Options::clear() {
   mPairs.clear();
}

std::string Options::toString() const {
   std::stringstream ss;
   for(int i = 0; i < mPairs.size(); i++) {
      if(i != 0)
         ss << " ";
      ss << mPairs[i].first << "=" << mPairs[i].second;
   }
   return ss.str();
}

Options::Options(std::string iOptionString) {
   addOptions(iOptionString);
}

void Options::addOption(std::string iOptionString) {
   int nextEqual = iOptionString.find('=');
   if(nextEqual < 0) {
      return;
   }
   else {
      std::string key   = iOptionString.substr(0, nextEqual);
      std::string value = iOptionString.substr(nextEqual+1);

      addOption(key, value);
   }
}
