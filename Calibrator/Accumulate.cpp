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
   std::vector<const Field*> inputFields(nTime);
   std::vector<Field*> outputFields(nTime);
   for(int t = 0; t < nTime; t++) {
      inputFields[t]  = &iFile.getField(mType, t);
      outputFields[t] = &iFile.getEmptyField();
   }

   // Parallelizable
   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         for(int e = 0; e < nEns; e++) {
            (*outputFields[0])[i][j][e] = 0;
            for(int t = 1; t < nTime; t++) {
               float previous = (*outputFields[t-1])[i][j][e];
               float current  = (*inputFields[t])[i][j][e];
               if(Util::isValid(current) && Util::isValid(previous)) {
                  (*outputFields[t])[i][j][e] = current + previous;
               }
               else {
                  (*outputFields[t])[i][j][e] = Util::MV;
               }
            }
         }
      }
   }

   // Add to file
   for(int t = 0; t < nTime; t++) {
      iFile.addField(*outputFields[t], mType, t);
   }
}
