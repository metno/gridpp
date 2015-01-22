#ifndef CALIBRATOR_H
#define CALIBRATOR_H
#include <string>
#include <vector>
class File;
class Options;

//! Abstract calibration class
class Calibrator {
   public:
      Calibrator();
      // Calibrate one or more fields
      void calibrate(File& iFile) const;
      static Calibrator* getScheme(std::string iType, Options& iOptions);

      // Ensure that ensemble members in iAfter are in the same order as in iBefore
      static void  shuffle(const std::vector<float>& iBefore, std::vector<float>& iAfter);
      virtual std::string name() const = 0;
   protected:
      virtual void calibrateCore(File& iFile) const = 0;
   private:
};
#include "Zaga.h"
#include "Cloud.h"
#include "Wind.h"
#include "Temperature.h"
#include "Accumulate.h"
#include "Smooth.h"
#endif
