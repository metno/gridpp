#include "Wind.h"
#include "../Util.h"

CalibratorWind::CalibratorWind(const ParameterFileRegion& iParameterFile) : Calibrator(iParameterFile) {

}
void CalibratorWind::calibrateCore(const DataFile& iInput, DataFile& iOutput) const {
   int nLat = iInput.getNumLat();
   int nLon = iInput.getNumLon();
   int nEns = iInput.getNumEns();
   int nTime = iInput.getNumTime();

   // Initialize calibrated fields in the output file
   std::vector<Field*> windCals(nTime);
   for(int t = 0; t < nTime; t++) {
      windCals[t] = &iOutput.getEmptyField();
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Parameters parameters = mParameterFile.getParameters(t);
      const Field& u = iInput.getField(Variable::U, t);
      const Field& v = iInput.getField(Variable::V, t);

      Field& windCal  = *windCals[t];

      // Parallelizable
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               float currU = u[i][j][e];
               float currV = v[i][j][e];
               if(Util::isValid(currU) && Util::isValid(currV)) {
                  float dir = getDir(currU, currV);
                  float speed = getSpeed(currU, currV);
                  float calSpeed = calibrate(speed, dir, parameters);
                  windCal[i][j][e] = calSpeed;
               }
            }
         }
      }
   }

   // Add to file
   for(int t = 0; t < nTime; t++) {
      iOutput.addField(*windCals[t], Variable::W, t);
   }
}

float CalibratorWind::calibrate(float iSpeed, float iDir, const Parameters& iParameters) const {

   return Util::MV;
}
