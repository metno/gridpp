#include "Wind.h"
#include "../Util.h"
#include "../File/File.h"

CalibratorWind::CalibratorWind(const ParameterFileRegion& iParameterFile) :
      Calibrator(),
      mParameterFile(iParameterFile) {
}
void CalibratorWind::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Initialize calibrated fields in the output file
   std::vector<FieldPtr> windCals(nTime);
   for(int t = 0; t < nTime; t++) {
      windCals[t] = iFile.getEmptyField();
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      FieldPtr u = iFile.getField(Variable::U, t);
      FieldPtr v = iFile.getField(Variable::V, t);

      Field& windCal  = *windCals[t];

      // Parallelizable
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            Parameters parameters = mParameterFile.getParameters(i, j, t);
            for(int e = 0; e < nEns; e++) {
               float currU = (*u)[i][j][e];
               float currV = (*v)[i][j][e];
               if(Util::isValid(currU) && Util::isValid(currV)) {
                  float dir = getDir(currU, currV);
                  float speed = getSpeed(currU, currV);
                  float calSpeed = calibrateLocal(speed, dir, parameters);
                  windCal[i][j][e] = calSpeed;
               }
            }
         }
      }
   }
}

float CalibratorWind::calibrateLocal(float iSpeed, float iDir, const Parameters& iParameters) const {
   return 16.4;
}

float CalibratorWind::getDir(float iU, float iV) {
   return 0;
}
float CalibratorWind::getSpeed(float iU, float iV) {
   return 0;
}
