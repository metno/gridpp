#ifndef CALIBRATION_H
#define CALIBRATION_H
#include "DataFile.h"
#include "ParameterFile.h"

class Calibration {
   public:
      Calibration(const ParameterFile& iParameterFile);
      void calibrate(const DataFile& iInput, DataFile& iOutput) const;
   private:
      float getCdf(float iX, float iEnsMean, float iEnsFrac) const;
      float getPdf(float iX, float iEnsMean, float iEnsFrac) const;
      float getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, Parameters& iParameters) const;
      float getP0(float iEnsMean, float iEnsFrac, Parameters& iParameters) const;
      const ParameterFile& mParameterFile;
      static const float mMaxEnsMean = 100;
      std::string mPrecipName;
      std::string mCloudName;
};
#endif
