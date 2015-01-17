#ifndef UTIL_H
#define UTIL_H
#include <string>
#include <vector>

typedef std::vector<std::vector<float> > vec2; // Lat, Lon
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
      static float getDistance(float lat1, float lon1, float lat2, float lon2);
      static float deg2rad(float deg);
      static float rad2deg(float rad);
      static float pi;
      static double radiusEarth;
   private:
      static bool mShowError;
      static bool mShowWarning;
      static bool mShowStatus;
};
#endif
