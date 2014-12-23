#ifndef CALIBRATION_H
#define CALIBRATION_H
#include "DataFile.h"
#include "ParameterFile.h"

//! Abstract calibration class
class Calibration {
   public:
      Calibration(const ParameterFile& iParameterFile);
      virtual void calibrate(const DataFile& iInput, DataFile& iOutput) const = 0;

      // Helper function
      static float logit(float p);
      static float invLogit(float x);
   protected:
      const ParameterFile& mParameterFile;
      template<class T1, class T2> struct sort_pair_second {
         bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.second < right.second;
         }
      };
      template<class T1, class T2> struct sort_pair_first {
         bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
         }
      };
};

class CalibrationPrecip : public Calibration {
   public:
      CalibrationPrecip(const ParameterFile& iParameterFile);
      void calibrate(const DataFile& iInput, DataFile& iOutput) const;
      //! Get probability mass at 0 mm (i.e probability of no precipitation)
      static float getP0(float iEnsMean, float iEnsFrac, Parameters& iParameters);
      //! Get Precipitation amount corresponding to quantile
      static float getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, int iTime, Parameters& iParameters);
      //! Set threshold between precip/no-precip used when computing what fraction of
      //! members have precip.
      void setFracThreshold(float iFraction);
   private:
      static const float mMaxEnsMean = 100;
      //! What precip threshold should be used to count members with no precip?
      float mFracThreshold;
      int mNeighbourhoodSize;
};
#endif
