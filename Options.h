#ifndef OPTIONS_H
#define OPTIONS_H
#include <iostream>
#include <vector>
#include <sstream>

class Options {
   public:
      Options();
      void addOption(std::string iKey, std::string iValue);
      template <class T> void addOption(std::string iKey, T iValue) {
         mKeys.push_back(iKey);
         std::stringstream ss;
         ss << iValue;
         mValues.push_back(ss.str());
      }
      void addOption(std::string iOptionString);
      void clear();
      template <class T> bool getValue(std::string iKey, T& iValue) {
         for(int i = 0; i < mKeys.size(); i++) {
            if(mKeys[i] == iKey) {
               std::stringstream ss (mValues[i]);
               ss >> iValue;
               return true;
            }
         }
         return false;
      };
   protected:
      std::vector<std::string> mKeys;
      std::vector<std::string> mValues;
};
#endif
