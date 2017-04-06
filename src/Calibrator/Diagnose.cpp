#include "Diagnose.h"
#include "../Util.h"
#include "../File/File.h"
CalibratorDiagnose::CalibratorDiagnose(Variable::Type iVariable, const Options& iOptions) :
      Calibrator(iOptions),
      mOutputVariable(iVariable) {
   /*
   std::vector<std::string> usingVariables;
   iOptions.getValues("using", usingVariables);
   for(int i = 0; i < usingVariables.size(); i++) {
      mDiagVariables.push_back(Variable::getType(usingVariables[i]));
   }
   */
}
bool CalibratorDiagnose::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Check that we have the required variables
   std::vector<Variable::Type> requiredVariables;
   if(mOutputVariable == Variable::W || mOutputVariable == Variable::WD) {
      if(iFile.hasVariable(Variable::U)) {
         requiredVariables.push_back(Variable::U);
         requiredVariables.push_back(Variable::V);
      }
      else {
         requiredVariables.push_back(Variable::Xwind);
         requiredVariables.push_back(Variable::Ywind);
      }
   }
   else if(mOutputVariable == Variable::U || mOutputVariable == Variable::V) {
      requiredVariables.push_back(Variable::W);
      requiredVariables.push_back(Variable::WD);
   }
   else if(mOutputVariable == Variable::Xwind || mOutputVariable == Variable::Ywind) {
      requiredVariables.push_back(Variable::W);
      requiredVariables.push_back(Variable::WD);
   }
   else {
      std::stringstream ss;
      ss << "Cannot diagnose " << Variable::getTypeName(mOutputVariable)
         << " because no conversion has been implemented";
      Util::error(ss.str());
   }

   for(int i = 0; i < requiredVariables.size(); i++) {
      if(!iFile.hasVariable(requiredVariables[i])) {
         std::stringstream ss;
         ss << "Cannot diagnose " << Variable::getTypeName(mOutputVariable)
            << " because " << Variable::getTypeName(requiredVariables[i]) << " is missing";
         Util::error(ss.str());
      }
   }

   // Get all fields
   for(int t = 0; t < nTime; t++) {
      Field& output = *iFile.getField(mOutputVariable, t);
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               // Diagnose wind speed from U and V
               if(mOutputVariable == Variable::W && iFile.hasVariable(Variable::U)) {
                  const Field& fieldU = *iFile.getField(Variable::U, t);
                  const Field& fieldV = *iFile.getField(Variable::V, t);
                  if(Util::isValid(fieldU(i,j,e)) && Util::isValid(fieldV(i,j,e)))
                     output(i,j,e) = sqrt(fieldU(i,j,e)*fieldU(i,j,e) + fieldV(i,j,e)*fieldV(i,j,e));
               }
               else if(mOutputVariable == Variable::W && iFile.hasVariable(Variable::Xwind)) {
                  const Field& fieldXwind = *iFile.getField(Variable::Xwind, t);
                  const Field& fieldYwind = *iFile.getField(Variable::Ywind, t);
                  if(Util::isValid(fieldXwind(i,j,e)) && Util::isValid(fieldYwind(i,j,e)))
                     output(i,j,e) = sqrt(fieldXwind(i,j,e)*fieldXwind(i,j,e) + fieldYwind(i,j,e)*fieldYwind(i,j,e));
               }
               // Diagnose U from W and WD
               else if(mOutputVariable == Variable::U || mOutputVariable == Variable::Xwind) {
                  const Field& fieldW = *iFile.getField(Variable::W, t);
                  const Field& fieldWD = *iFile.getField(Variable::WD, t);
                  if(Util::isValid(fieldW(i,j,e)) && Util::isValid(fieldWD(i,j,e)))
                     output(i,j,e) = - fieldW(i,j,e) * sin(fieldWD(i,j,e) / 180.0 * Util::pi);
               }
               // Diagnose V from W and WD
               else if(mOutputVariable == Variable::V || mOutputVariable == Variable::Ywind) {
                  const Field& fieldW = *iFile.getField(Variable::W, t);
                  const Field& fieldWD = *iFile.getField(Variable::WD, t);
                  if(Util::isValid(fieldW(i,j,e)) && Util::isValid(fieldWD(i,j,e)))
                     output(i,j,e) = - fieldW(i,j,e) * cos(fieldWD(i,j,e) / 180.0 * Util::pi);
               }
            }
         }
      }
   }
   return true;
}
std::string CalibratorDiagnose::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c diagnose","Will attempt to diagnose the variable by inspecting the other available variables. Currently only W, WD, U, and V are supported. This is useful when one variable has been post-processed leaving other variables that are related to these in need of update.") << std::endl;
   return ss.str();
}
