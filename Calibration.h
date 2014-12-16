#ifndef CALIBRATION_H
#define CALIBRATION_H
#include "DataFile.h"
#include "ParameterFile.h"

class Calibration {
   public:
      Calibration(const ParameterFile& iParameterFile);
      void calibrate(const DataFile& iInput, DataFile& iOutput) const;
      static float logit(float p);
      static float invLogit(float x);
      static float getP0(float iEnsMean, float iEnsFrac, Parameters& iParameters);
      static float getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, int iTime, Parameters& iParameters);
   private:
      float getCdf(float iX, float iEnsMean, float iEnsFrac) const;
      float getPdf(float iX, float iEnsMean, float iEnsFrac) const;
      const ParameterFile& mParameterFile;
      static const float mMaxEnsMean = 100;
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
      //! What precip threshold should be used to count members with no precip?
      float mFracThreshold;
};
#endif
