#ifndef CALIBRATOR_H
#define CALIBRATOR_H
#include "../DataFile.h"
#include "../ParameterFile.h"

//! Abstract calibration class
class Calibrator {
   public:
      Calibrator(const ParameterFile& iParameterFile);
      // Calibrate one or more fields
      void calibrate(DataFile& iFile) const;

      // Helper function
      static float logit(float p);
      static float invLogit(float x);
      // Ensure that ensemble members in iAfter are in the same order as in iBefore
      static void  shuffle(const std::vector<float>& iBefore, std::vector<float>& iAfter);
   protected:
      virtual void calibrateCore(DataFile& iFile) const = 0;
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
   private:
};
#include "Zaga.h"
#include "Cloud.h"
//#include "Wind.h"
#include "Temperature.h"
#include "Accumulate.h"
#endif
