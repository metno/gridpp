#ifndef UTIL_H
#define UTIL_H
#include <string>

class Util {
   public:
      Util();
      //! Aborts the program with an error message
      static void error(std::string iMessage);
      static void warning(std::string iMessage);
      static void status(std::string iMessage);

      //! Returns the current unix time in seconds
      static double clock();
      static void setShowError(bool flag);
      static void setShowWarning(bool flag);
      static void setShowStatus(bool flag);
      //! Missing value indicator
      static float MV;
      //! Checks if a value is valid
      static bool isValid(float iValue);
   private:
      static bool mShowError;
      static bool mShowWarning;
      static bool mShowStatus;
};
#endif
