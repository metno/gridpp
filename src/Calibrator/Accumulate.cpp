#include "Accumulate.h"
#include "../Util.h"
#include "../File/File.h"
CalibratorAccumulate::CalibratorAccumulate(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mOutputVariable(Variable::PrecipAcc) {
   std::string outputVariable;
   if(iOptions.getValue("outputVariable", outputVariable)) {
      mOutputVariable = Variable::getType(outputVariable);
   }
}
bool CalibratorAccumulate::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   if(!iFile.hasVariable(mVariable)) {
      Util::error("File '" + iFile.getFilename() + "' does not have variable '" + mVariable.getName() + "'-");
   }

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

      Variable variable;
      bool found = iFile.getVariable(mOutputVariable, variable);
      if(!found)
         Util::error("File does not have accumulated field defined");
      iFile.addField(fieldsAcc[t], variable, t);
   }
   return true;
}
std::string CalibratorAccumulate::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c accumulate","Accumlates a value over time. Used to accumulate precipitation. It is assumed that the raw value for time t is the precip on the interval [t-1,t]. After accumulation, the value for time t is then the sum over the interval [0,t]. Thus, the accumulated value for the first timestep will be missing and the raw value for time t=0 will never be used.") << std::endl;
   ss << Util::formatDescription("   outputVariable=PrecipAcc","Which variable should the accumulate be written to?") << std::endl;
   return ss.str();
}
