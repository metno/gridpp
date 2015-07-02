#include "Window.h"
#include "../Util.h"
#include "../File/File.h"
CalibratorWindow::CalibratorWindow(Variable::Type iVariable, const Options& iOptions) :
      Calibrator(NULL, iOptions),
      mRadius(3),
      mVariable(iVariable),
      mStatType(Util::StatTypeMean),
      mQuantile(Util::MV) {
   iOptions.getValue("radius", mRadius);
   if(mRadius < 0) {
      std::stringstream ss;
      ss << "CalibratorWindow: 'radius' (" << mRadius << ") must be >= 0";
      Util::error(ss.str());
   }

   std::string op;
   if(iOptions.getValue("stat", op)) {
      if(op == "mean") {
         mStatType = Util::StatTypeMean;
      }
      else if(op == "min") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 0;
      }
      else if(op == "max") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 1;
      }
      else if(op == "median") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 0.5;
      }
      else if(op == "std") {
         mStatType = Util::StatTypeStd;
      }
      else if(op == "quantile"){
         mStatType = Util::StatTypeQuantile;
         if(!iOptions.getValue("quantile", mQuantile)) {
            Util::error("CalibratorWindow: option 'quantile' is required");
         }
         if(!Util::isValid(mQuantile) || mQuantile < 0 || mQuantile > 1) {
            Util::error("CalibratorWindow: 'quantile' must be on the interval [0,1]");
         }
      }
      else {
         Util::error("CalibratorWindow: Unrecognized value for 'stat'");
      }
   }
}
bool CalibratorWindow::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Get all fields
   std::vector<FieldPtr> fields(nTime);
   std::vector<FieldPtr> fieldsOrig(nTime);
   for(int t = 0; t < nTime; t++) {
      fieldsOrig[t]     = iFile.getField(mVariable, t);
      fields[t] = iFile.getEmptyField();
   }

   for(int t = 0; t < nTime; t++) {
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               std::vector<float> window;
               window.resize(2*mRadius+1, Util::MV);
               int index = 0;
               for(int tt = std::max(0,t-mRadius); tt <= std::min(nTime-1, t + mRadius); tt++) {
                  float curr = (*fieldsOrig[tt])(i,j,e);
                  window[index] = curr;
                  index++;
               }
               (*fields[t])(i,j,e) = Util::calculateStat(window, mStatType, mQuantile);
            }
         }
      }
   }
   for(int t = 0; t < nTime; t++) {
      iFile.addField(fields[t], mVariable, t);
   }
   return true;
}

std::string CalibratorWindow::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c window","Applies a statistical operator to values within a temporal window. Any missing values are ignored when computing the statistic.") << std::endl;
   ss << Util::formatDescription("   radius=required", "Define the window as all offsets within +- radius (must be 0 or greater). The unit is in number of time indices in the file.") << std::endl;
   ss << Util::formatDescription("   stat=mean", "What statistical operator should be applied to the window? One of 'mean', 'median', 'min', 'max', 'std', or 'quantile'. 'std' is the population standard deviation.") << std::endl;
   ss << Util::formatDescription("   quantile=undef", "If stat=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
   return ss.str();
}
