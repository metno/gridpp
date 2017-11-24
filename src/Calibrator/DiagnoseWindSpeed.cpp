#include "DiagnoseWindSpeed.h"
#include "../Util.h"
#include "../File/File.h"
#include <cmath>
CalibratorDiagnoseWindSpeed::CalibratorDiagnoseWindSpeed(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
   iOptions.getRequiredValue("x", mXVariable);
   iOptions.getRequiredValue("y", mYVariable);
}
bool CalibratorDiagnoseWindSpeed::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Get all fields
   for(int t = 0; t < nTime; t++) {
      Field& output = *iFile.getField(mVariable, t);
      const Field& fieldX = *iFile.getField(mXVariable, t);
      const Field& fieldY = *iFile.getField(mYVariable, t);

      #pragma omp parallel for
      for(int y = 0; y < nY; y++) {
         for(int x = 0; x < nX; x++) {
            for(int e = 0; e < nEns; e++) {
               if(Util::isValid(fieldX(y, x, e)) && Util::isValid(fieldY(y, x, e)))
                  output(y, x, e) = sqrt(fieldX(y, x, e)*fieldX(y, x, e) + fieldY(y, x, e)*fieldY(y, x, e));
            }
         }
      }
   }
   return true;
}
std::string CalibratorDiagnoseWindSpeed::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c diagnoseWindSpeed","Compute wind speed from x, y components.") << std::endl;
   ss << Util::formatDescription("   x","X-wind variable name") << std::endl;
   ss << Util::formatDescription("   y","Y-wind variable name") << std::endl;
   return ss.str();
}

