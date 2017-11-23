#include "Diagnose.h"
#include "../Util.h"
#include "../File/File.h"
#include <cmath>
CalibratorDiagnose::CalibratorDiagnose(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
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
   if(mVariable.getType() == Variable::W || mVariable.getType() == Variable::WD) {
      if(iFile.hasVariable(Variable::U)) {
         requiredVariables.push_back(Variable::U);
         requiredVariables.push_back(Variable::V);
      }
      else {
         requiredVariables.push_back(Variable::Xwind);
         requiredVariables.push_back(Variable::Ywind);
      }
   }
   else if(mVariable.getType() == Variable::U || mVariable.getType() == Variable::V) {
      requiredVariables.push_back(Variable::W);
      requiredVariables.push_back(Variable::WD);
   }
   else if(mVariable.getType() == Variable::Xwind || mVariable.getType() == Variable::Ywind) {
      requiredVariables.push_back(Variable::W);
      requiredVariables.push_back(Variable::WD);
   }
   else if(mVariable.getType() == Variable::RH) {
      requiredVariables.push_back(Variable::T);
      requiredVariables.push_back(Variable::TD);
   }
   else if(mVariable.getType() == Variable::TD) {
      requiredVariables.push_back(Variable::T);
      requiredVariables.push_back(Variable::RH);
   }
   else {
      std::stringstream ss;
      ss << "Cannot diagnose " << Variable::getTypeName(mVariable.getType())
         << " because no conversion has been implemented";
      Util::error(ss.str());
   }

   for(int i = 0; i < requiredVariables.size(); i++) {
      if(!iFile.hasVariableWithoutDeriving(requiredVariables[i])) {
         std::stringstream ss;
         ss << "Cannot diagnose " << Variable::getTypeName(mVariable.getType())
            << " because " << Variable::getTypeName(requiredVariables[i]) << " is missing";
         Util::error(ss.str());
      }
   }

   // Get all fields
   for(int t = 0; t < nTime; t++) {
      Field& output = *iFile.getField(mVariable, t);
      // Diagnose W from U and V
      if(mVariable.getType() == Variable::W && iFile.hasVariableWithoutDeriving(Variable::U)) {
         const Field& fieldU = *iFile.getField(Variable::U, t);
         const Field& fieldV = *iFile.getField(Variable::V, t);
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(fieldU(i,j,e)) && Util::isValid(fieldV(i,j,e)))
                     output(i,j,e) = sqrt(fieldU(i,j,e)*fieldU(i,j,e) + fieldV(i,j,e)*fieldV(i,j,e));
               }
            }
         }
      }
      // Diagnose W from Xwind and Ywind
      else if(mVariable.getType() == Variable::W && iFile.hasVariableWithoutDeriving(Variable::Xwind)) {
         const Field& fieldXwind = *iFile.getField(Variable::Xwind, t);
         const Field& fieldYwind = *iFile.getField(Variable::Ywind, t);
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(fieldXwind(i,j,e)) && Util::isValid(fieldYwind(i,j,e))) {
                      output(i,j,e) = sqrt(fieldXwind(i,j,e)*fieldXwind(i,j,e) + fieldYwind(i,j,e)*fieldYwind(i,j,e));
                  }
               }
            }
         }
      }
      // Diagnose WD from U and V
      else if(mVariable.getType() == Variable::WD && iFile.hasVariableWithoutDeriving(Variable::U)) {
         const Field& fieldU = *iFile.getField(Variable::U, t);
         const Field& fieldV = *iFile.getField(Variable::V, t);
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(fieldU(i,j,e)) && Util::isValid(fieldV(i,j,e))) {
                     float dir = std::atan2(-fieldU(i,j,e), -fieldV(i,j,e)) * 180 / Util::pi;
                     if(dir < 0)
                        dir += 360;
                     output(i,j,e) = dir;
                  }
               }
            }
         }
      }
      // Diagnose WD from Xwind and Ywind
      else if(mVariable.getType() == Variable::WD && iFile.hasVariableWithoutDeriving(Variable::Xwind)) {
         const Field& fieldXwind = *iFile.getField(Variable::Xwind, t);
         const Field& fieldYwind = *iFile.getField(Variable::Ywind, t);
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(fieldXwind(i,j,e)) && Util::isValid(fieldYwind(i,j,e))) {
                     float dir = std::atan2(-fieldXwind(i,j,e), -fieldYwind(i,j,e)) * 180 / Util::pi;
                     if(dir < 0)
                        dir += 360;
                     output(i,j,e) = dir;
                  }
               }
            }
         }
      }
      // Diagnose U from W and WD
      else if(mVariable.getType() == Variable::U || mVariable.getType() == Variable::Xwind) {
         const Field& fieldW = *iFile.getField(Variable::W, t);
         const Field& fieldWD = *iFile.getField(Variable::WD, t);
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(fieldW(i,j,e)) && Util::isValid(fieldWD(i,j,e))) {
                     output(i,j,e) = -fieldW(i,j,e) * sin(fieldWD(i,j,e) / 180.0 * Util::pi);
                  }
               }
            }
         }
      }
      // Diagnose V from W and WD
      else if(mVariable.getType() == Variable::V || mVariable.getType() == Variable::Ywind) {
         const Field& fieldW = *iFile.getField(Variable::W, t);
         const Field& fieldWD = *iFile.getField(Variable::WD, t);
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(fieldW(i,j,e)) && Util::isValid(fieldWD(i,j,e))) {
                     output(i,j,e) = - fieldW(i,j,e) * cos(fieldWD(i,j,e) / 180.0 * Util::pi);
                  }
               }
            }
         }
      }
      // Diagnose RH from T and TD
      else if(mVariable.getType() == Variable::RH) {
         const Field& fieldT = *iFile.getField(Variable::T, t);
         const Field& fieldTD = *iFile.getField(Variable::TD, t);
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(fieldT(i,j,e)) && Util::isValid(fieldTD(i,j,e))) {
                     output(i,j,e) = dewpoint2RH(fieldT(i,j,e), fieldTD(i,j,e));
                  }
               }
            }
         }
      }
      // Diagnose TD from T and RH
      else if(mVariable.getType() == Variable::TD) {
         const Field& fieldT = *iFile.getField(Variable::T, t);
         const Field& fieldRH = *iFile.getField(Variable::RH, t);
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(fieldT(i,j,e)) && Util::isValid(fieldRH(i,j,e))) {
                     output(i,j,e) = RH2dewpoint(fieldT(i,j,e), fieldRH(i,j,e));
                  }
               }
            }
         }
      }
      else {
         // This part should never happen, since it should have been caught ealier
         std::stringstream ss;
         ss << "Cannot diagnose " << Variable::getTypeName(mVariable.getType());
         Util::error(ss.str());
      }
   }
   return true;
}

float CalibratorDiagnose::mEwt[41] = { .000034,.000089,.000220,.000517,.001155,.002472,
                            .005080,.01005, .01921, .03553, .06356, .1111,
                            .1891,  .3139,  .5088,  .8070,  1.2540, 1.9118,
                            2.8627, 4.2148, 6.1078, 8.7192, 12.272, 17.044,
                            23.373, 31.671, 42.430, 56.236, 73.777, 95.855,
                            123.40, 157.46, 199.26, 250.16, 311.69, 385.56,
                            473.67, 578.09, 701.13, 845.28, 1013.25 };

float CalibratorDiagnose::dewpoint2RH(float iTemperature, float iDewPointTemperature) {
   // Taken from https://github.com/metno/wdb2ts
   if(Util::isValid(iTemperature) && Util::isValid(iDewPointTemperature)) {
      float x, et, etd, rh;
      int l;

      x = (iTemperature - 173.16) * 0.2;
      l = int(x);
      et = mEwt[l] + (mEwt[l + 1] - mEwt[l]) * (x - float(l));

      x = (iDewPointTemperature - 173.16) * 0.2;
      l = int(x);
      etd = mEwt[l] + (mEwt[l + 1] - mEwt[l]) * (x - float(l));
      rh = etd / et;

      if(rh < 0)
         rh = 0;
      if(rh > 1)
         rh = 1;

      return rh;
   }
   else
      return Util::MV;
}


float CalibratorDiagnose::RH2dewpoint(float iTemperature, float iRelativeHumidity) {
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


std::string CalibratorDiagnose::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c diagnose","Will attempt to diagnose the variable by inspecting the other available variables. Currently W, WD, U, V, Xwind, Ywind, RH, and TD are supported. This is useful when one variable has been post-processed leaving other variables that are related to these in need of update.") << std::endl;
   return ss.str();
}

