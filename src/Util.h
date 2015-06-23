#ifndef UTIL_H
#define UTIL_H
#include <string>
#include <vector>
#include <set>

typedef std::vector<std::vector<float> > vec2; // Lat, Lon

//! Contains helper functions
class Util {
   public:
      //! Aborts the program with an error message to std::cout
      static void error(std::string iMessage);

      //! Issues a warning message to std::cout
      static void warning(std::string iMessage);

      //! Issues a status message to std::cout
      static void status(std::string iMessage);

      //! True will force calls to Util::error to show error messages.
      //! Calls to errors will cause abort regardless of this setting.
      static void setShowError(bool flag);
      
      //! True will force calls to Util::warning to show warning messages
      static void setShowWarning(bool flag);
      
      //! True will force calls to Util::status to show status messages
      static void setShowStatus(bool flag);
      
      //! Missing value indicator
      static float MV;

      //! \brief Checks if a value is valid
      //! @return false if iValue equals Util::MV or is +-inf or nan. Returns true otherwise.
      static bool isValid(float iValue);

      //! Checks if the file exists
      static bool exists(const std::string& iFilename);

      //! \brief Computes the geographical distance between two locations, assuming a spherical earth.
      //! @param lat1 Latitude in degrees of point 1
      //! @param lon1 Longitude in degrees of point 1
      //! @param lat2 Latitude in degrees of point 2
      //! @param lon2 Longitude in degrees of point 2
      //! @param approx Should an approximation be used? If true, assumes a flat earth in the
      //!               vicinity of the two points. Computes the result much faster.
      //! @return Distance in meters
      static float getDistance(float lat1, float lon1, float lat2, float lon2, bool approx=false);

      //! Returns the current unix time in seconds
      static double clock();
      
      //! Convert degrees to radians
      static float deg2rad(float deg);
      
      //! Convert radians to degrees
      static float rad2deg(float rad);

      static float pi;

      static std::string gridppVersion();
      
      //! Radius of the earth in meters
      static double radiusEarth;
      
      //! \brief Convert unix timestamp to YYYYMMDD format
      static int getDate(time_t iUnixTime);

      //! \brief Convert unix timestamp to HHMMSS format
      static int getTime(time_t iUnixTime);
      
      //! \brief Get the current unix time (UTC)
      //! @param iDate in YYYYMMDD format
      //! @param iTime in HHMMSS format (can be negative to subtract time)
      static time_t getUnixTime(int iDate, int iTime);

      //! \brief Get the current date in YYYYMMDD format (UTC)
      static int getCurrentDate();

      //! \brief Get the current date in HHMMSS format (UTC)
      static int getCurrentTime();

      //! \brief Get the current date/time in unix time (UTC)
      static time_t getCurrentUnixTime();

      //! \brief Get the current date in YYYY-MM-DD HH:MM:SS format
      static std::string getCurrentTimeStamp();
      
      //! \brief Calculates new date by adding hours to an existing date
      //! @param iDate Date in format YYYYMMDD
      //! @param iAddHours Number of hours to add
      //! @return New date in format YYYYMMDD
      static int calcDate(int iDate, float iAddHours);
      
      //! \brief Creates a vector of strings by splitting string on each space
      //! Multiple consecutive spaces are treated as one space
      static std::vector<std::string> split(std::string iString);

      //! \brief Computes the logit of p
      //! @param p Must lie in the interval (0,1)
      static float logit(float p);

      //! \brief Computes the inverse logit of x
      //! @return x A value on the interval (0,1)
      static float invLogit(float x);

      template <class T> static std::vector<T> combine(const std::vector<T>& i1, const std::vector<T>& i2) {
         std::set<T> allValues(i1.begin(), i1.end());
         for(int i = 0; i < i2.size(); i++) {
            allValues.insert(i2[i]);
         }
         return std::vector<T>(allValues.begin(), allValues.end());
      };
      //! Returns true if iString contains one or more of iChar
      static bool hasChar(std::string iString, char iChar);

      //! Create a formatted string useful for outputing descriptions of methods. Places the title
      //! on the left side and the message on the right, with automatic line wrapping.
      static std::string formatDescription(std::string iTitle, std::string iMessage, int iTitleLength=32, int iMaxLength=100, int iTitleIndent=3);

      enum OperatorType {
         OperatorMean      = 0,
         OperatorStd       = 30,
         OperatorQuantile  = 40
      };

      //! Applies statistics operator to array. Missing values are ignored.
      static float applyOperator(const std::vector<float>& iArray, Util::OperatorType iOperator, float iQuantile=Util::MV);
      
      //! \brief Comparator class for sorting pairs using the first entry.
      //! Sorts from smallest to largest
      template<class T1, class T2> struct sort_pair_first {
         bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
         }
      };

      //! Computes the matrix inverse
      //! If any elements in iMatrix is missing, then all elements in the inverse will be missing.
      //! If the matrix is not invertable, then all elements in the inverse will be missing.
      static vec2 inverse(const vec2 iMatrix);

      //! Copy the file with filename iFrom to filename iTo
      static bool copy(std::string iFrom, std::string iTo);
     
      //! \brief Comparator class for sorting pairs using the second entry.
      //! Sorts from smallest to largest
      template<class T1, class T2> struct sort_pair_second {
         bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.second < right.second;
         }
      };
   private:
      static bool mShowError;
      static bool mShowWarning;
      static bool mShowStatus;
};
#endif
