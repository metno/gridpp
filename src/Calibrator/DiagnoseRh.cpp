#include "DiagnoseRh.h"
#include "../Util.h"
#include "../File/File.h"
#include <cmath>
CalibratorDiagnoseRh::CalibratorDiagnoseRh(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
   iOptions.getRequiredValue("temperature", mTemperature);
   iOptions.getRequiredValue("dewpoint", mDewpoint);
   iOptions.check();
}
bool CalibratorDiagnoseRh::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Get all fields
   for(int t = 0; t < nTime; t++) {
      Field& output = *iFile.getField(mVariable, t);
      const Field& fieldT = *iFile.getField(mTemperature, t);
      const Field& fieldTd = *iFile.getField(mDewpoint, t);

      #pragma omp parallel for
      for(int y = 0; y < nY; y++) {
         for(int x = 0; x < nX; x++) {
            for(int e = 0; e < nEns; e++) {
               if(Util::isValid(fieldT(y, x ,e)) && Util::isValid(fieldTd(y, x, e))) {
                  output(y, x, e) = dewpoint2RH(fieldT(y, x, e), fieldTd(y, x, e));
               }
            }
         }
      }
   }
   return true;
}

float CalibratorDiagnoseRh::mEwt[41] = { .000034,.000089,.000220,.000517,.001155,.002472,
                            .005080,.01005, .01921, .03553, .06356, .1111,
                            .1891,  .3139,  .5088,  .8070,  1.2540, 1.9118,
                            2.8627, 4.2148, 6.1078, 8.7192, 12.272, 17.044,
                            23.373, 31.671, 42.430, 56.236, 73.777, 95.855,
                            123.40, 157.46, 199.26, 250.16, 311.69, 385.56,
                            473.67, 578.09, 701.13, 845.28, 1013.25 };

float CalibratorDiagnoseRh::dewpoint2RH(float iTemperature, float iDewPointTemperature) {
   // Taken from https://github.com/metno/wdb2ts
   if(Util::isValid(iTemperature) && Util::isValid(iDewPointTemperature)) {
      float x, et, etd, rh;
      int l;

      x = (iTemperature - 173.16) * 0.2;
      l = int(x);
      et = mEwt[l] + (mEwt[l + 1] - mEwt[l]) * (x - float(l));

      x = (iDewPointTemperature - 173.16) * 0.2;
      l = int(x);
      etd = mEwt[l] + (mEwt[l + 1] - mEwt[l]) * (x - float(l));
      rh = etd / et;

      if(rh < 0)
         rh = 0;
      if(rh > 1)
         rh = 1;

      return rh;
   }
   else
      return Util::MV;
}



std::string CalibratorDiagnoseRh::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c diagnoseRh","Diagnoses dewpoint temperature.") << std::endl;
   ss << Util::formatDescription("   temperature","Name of temperature variable") << std::endl;
   ss << Util::formatDescription("   dewpoint","Name of dewpoint temperature variable") << std::endl;
   return ss.str();
}

