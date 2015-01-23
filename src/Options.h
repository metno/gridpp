#ifndef OPTIONS_H
#define OPTIONS_H
#include <iostream>
#include <vector>
#include <sstream>

// Container class for multiple key, value pairs
// Empty (i.e. "") keys or values are not allowed
class Options {
   public:
      Options(std::string iOptionString="");
      // Neither iKey nor iValue can be ""
      void addOption(std::string iKey, std::string iValue);
      // iOptionString is of the format: "key=value"
      // If key= is passed, nothing is added
      void addOption(std::string iOptionString);
      // Convert iValue to string and add to options
      template <class T> void addOption(std::string iKey, T iValue) {
         std::stringstream ss;
         ss << iValue;
         std::string value = ss.str();
         addOption(iKey, value);
      }
      // Remove all options
      void clear();
      // Returns true if key is found and stores in iValue. If not found,
      // iValue is unchanged.
      template <class T> bool getValue(std::string iKey, T& iValue) const {
         for(int i = 0; i < mKeys.size(); i++) {
            if(mKeys[i] == iKey) {
               std::stringstream ss (mValues[i]);
               ss >> iValue;
               return true;
            }
         }
         return false;
      };
      /* Why do I need this?
      bool getValue(std::string iKey, std::string& iValue) {
         for(int i = 0; i < mKeys.size(); i++) {
            if(mKeys[i] == iKey) {
               iValue = mValues[i];
               return true;
            }
         }
         return false;
      };
      */
   protected:
      std::vector<std::string> mKeys;
      std::vector<std::string> mValues;
};
#endif
