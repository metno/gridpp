#include "Accumulate.h"
#include "../Util.h"
#include "../File/File.h"
CalibratorAccumulate::CalibratorAccumulate(Variable::Type iVariable) :
      Calibrator(),
      mVariable(iVariable) {
}
bool CalibratorAccumulate::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Get all fields
   std::vector<FieldPtr> fields(nTime);
   std::vector<FieldPtr> fieldsAcc(nTime);
   for(int t = 0; t < nTime; t++) {
      fields[t]    = iFile.getField(mVariable, t);
      fieldsAcc[t] = iFile.getEmptyField();
   }

   for(int t = 0; t < nTime; t++) {
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               if(t == 0) {
                  (*fieldsAcc[t])(i,j,e) = 0;
               }
               else {
                  float previous = (*fieldsAcc[t-1])(i,j,e);
                  float current  = (*fields[t])(i,j,e);
                  if(Util::isValid(current) && Util::isValid(previous)) {
                     (*fieldsAcc[t])(i,j,e) = current + previous;
                  }
                  else {
                     (*fieldsAcc[t])(i,j,e) = Util::MV;
                  }
               }
            }
         }
      }

      Variable::Type variableAcc;
      if(mVariable == Variable::Precip) {
         variableAcc = Variable::PrecipAcc;
      }
      else {
         Util::error("Cannot accumulate " + Variable::getTypeName(mVariable));
      }
      iFile.addField(fieldsAcc[t], variableAcc, t);
   }
   return true;
}
std::string CalibratorAccumulate::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c accumulate","") << std::endl;
   return ss.str();
}
