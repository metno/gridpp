#ifndef OPTIONS_H
#define OPTIONS_H
#include <iostream>
#include <vector>
#include <sstream>

//! Container class for key-value pairs.
//! Empty (i.e. "") keys or values are not allowed
class Options {
   public:
      //! Creates container
      // @param iOptionString options with format: "key1=value1 key2=value2..."
      Options(std::string iOptionString="");

      //! \brief Adds key and value to container
      //! @param iKey cannot be ""
      //! @param iValue cannot be ""
      void addOption(std::string iKey, std::string iValue);

      //! \brief Parses options and adds to container
      //! @param iOptionString with format: "key1=value1 key2=value2...". key= ignored.
      void addOptions(std::string iOptionString);

      //! Convert iValue to string and add to options
      template <class T> void addOption(std::string iKey, T iValue) {
         std::stringstream ss;
         ss << iValue;
         std::string value = ss.str();
         addOption(iKey, value);
      }

      //! Remove all key-value pairs from container
      void clear();

      //! \brief Find value corresponding to key
      //! @param iKey find this key
      //! @param iValue places the value in this variable. If key does not exist, this is unchanged
      //! @return true if key is found, false otherwise
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
      //! Parse a string with a single option "key=value"
      void addOption(std::string iOptionString);
      std::vector<std::string> mKeys;
      std::vector<std::string> mValues;
};
#endif
