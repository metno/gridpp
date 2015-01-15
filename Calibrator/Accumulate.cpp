#include "Accumulate.h"
#include "../Util.h"
CalibratorAccumulate::CalibratorAccumulate(const ParameterFile& iParameterFile, Variable::Type iType) :
      Calibrator(iParameterFile),
      mType(iType) {
}
void CalibratorAccumulate::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Get all fields
   std::vector<FieldPtr> fields(nTime);
   std::vector<FieldPtr> fieldsAcc(nTime);
   for(int t = 0; t < nTime; t++) {
      fields[t]    = iFile.getField(mType, t);
      fieldsAcc[t] = iFile.getEmptyField();
   }

   for(int t = 0; t < nTime; t++) {
      // Parallelizable
#pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               if(t == 0) {
                  (*fieldsAcc[t])[i][j][e] = 0;
               }
               else {
                  float previous = (*fieldsAcc[t-1])[i][j][e];
                  float current  = (*fields[t])[i][j][e];
                  if(Util::isValid(current) && Util::isValid(previous)) {
                     (*fieldsAcc[t])[i][j][e] = current + previous;
                  }
                  else {
                     (*fieldsAcc[t])[i][j][e] = Util::MV;
                  }
               }
            }
         }
      }
      iFile.addField(fieldsAcc[t], Variable::PrecipAcc, t);
   }
}
