#include "DiagnoseDewpoint.h"
#include "../Util.h"
#include "../File/File.h"
#include <cmath>
CalibratorDiagnoseDewpoint::CalibratorDiagnoseDewpoint(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
   iOptions.getRequiredValue("temperature", mTemperature);
   iOptions.getRequiredValue("rh", mRH);
   iOptions.check();
}
bool CalibratorDiagnoseDewpoint::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Get all fields
   for(int t = 0; t < nTime; t++) {
      Field& output = *iFile.getField(mVariable, t);
      const Field& fieldT = *iFile.getField(mTemperature, t);
      const Field& fieldRH = *iFile.getField(mRH, t);

      #pragma omp parallel for
      for(int y = 0; y < nY; y++) {
         for(int x = 0; x < nX; x++) {
            for(int e = 0; e < nEns; e++) {
               if(Util::isValid(fieldT(y, x ,e)) && Util::isValid(fieldRH(y, x, e))) {
                  output(y, x, e) = RH2dewpoint(fieldT(y, x, e), fieldRH(y, x, e));
               }
            }
         }
      }
   }
   return true;
}

float CalibratorDiagnoseDewpoint::RH2dewpoint(float iTemperature, float iRelativeHumidity) {
   if(Util::isValid(iTemperature) && Util::isValid(iRelativeHumidity)) {
      // Taken from https://github.com/metno/wdb2ts
      float tempC = iTemperature - 273.15;
      float e = (iRelativeHumidity)*0.611*exp( (17.63 * tempC) / (tempC + 243.04) );
      float tdC = (116.9 + 243.04 * log( e ))/(16.78-log( e ));
      float td = tdC + 273.15;
      return (td<=iTemperature ? td : iTemperature);
      // Taken from https://github.com/WFRT/comps
      // float es = 0.611*exp(5423.0*(1/273.15 - 1/(iTemperature)));
      // float e  = iRelativeHumidity*(es);
      // float td = 1/(1/273.15 - 1.844e-4*log(e/0.611));
      // return td;
   }
   else
      return Util::MV;
}


std::string CalibratorDiagnoseDewpoint::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c diagnoseDewpoint","Diagnoses dewpoint temperature.") << std::endl;
   ss << Util::formatDescription("   temperature","Name of temperature variable") << std::endl;
   ss << Util::formatDescription("   rh","Name of relative humidity variable") << std::endl;
   return ss.str();
}

