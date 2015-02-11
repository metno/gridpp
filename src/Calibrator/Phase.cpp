#include "Phase.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorPhase::CalibratorPhase(const ParameterFile* iParameterFile) :
      Calibrator(),
      mParameterFile(iParameterFile),
      mMinPrecip(0.2),
      mEstimatePressure(true),
      mUseWetbulb(1) {
   if(iParameterFile->getNumParameters() != 2) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' does not have two datacolumns");
   }

}
bool CalibratorPhase::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   iFile.initNewVariable(Variable::Phase);
   vec2 elevs = iFile.getElevs();


   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      const Parameters& par = mParameterFile->getParameters(t);
      float snowSleetThreshold = par[0];
      float sleetRainThreshold = par[1];
      const FieldPtr temp = iFile.getField(Variable::T, t);
      const FieldPtr precip = iFile.getField(Variable::Precip, t);
      FieldPtr phase = iFile.getField(Variable::Phase, t);
      FieldPtr pressure;
      FieldPtr rh;
      if(mUseWetbulb) {
         // Only load these fields if they are to be used, to save memory
         rh = iFile.getField(Variable::RH, t);
         if(!mEstimatePressure)
            pressure = iFile.getField(Variable::P, t);
      }

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            float currElev = elevs[i][j];
            for(int e = 0; e < nEns; e++) {
               float currDryTemp  = (*temp)(i,j,e);
               float currTemp     = currDryTemp;
               float currPrecip   = (*precip)(i,j,e);
               if(mUseWetbulb) {
                  float currPressure;
                  if(mEstimatePressure)
                     currPressure = DownscalerPressure::calcPressure(0, 101325, currElev);
                  else
                     currPressure = (*pressure)(i,j,e);
                  float currRh       = (*rh)(i,j,e);
                  float currWetbulb  = getWetbulb(currDryTemp, currPressure, currRh);
                  if(Util::isValid(snowSleetThreshold) && Util::isValid(sleetRainThreshold)
                        && Util::isValid(currPrecip)   && Util::isValid(currDryTemp)
                        && Util::isValid(currPressure) && Util::isValid(currRh)) {
                     currTemp = currWetbulb;
                  }
               }
               if(Util::isValid(snowSleetThreshold) && Util::isValid(sleetRainThreshold) && Util::isValid(currTemp)) {
                  if(currPrecip <= mMinPrecip)
                     (*phase)(i,j,e)  = Variable::PhaseNone;
                  else if(!Util::isValid(currTemp))
                     (*phase)(i,j,e)  = Util::MV;
                  else if(currTemp <= snowSleetThreshold)
                     (*phase)(i,j,e)  = Variable::PhaseSnow;
                  else if(currTemp <= sleetRainThreshold)
                     (*phase)(i,j,e)  = Variable::PhaseSleet;
                  else
                     (*phase)(i,j,e)  = Variable::PhaseRain;
               }
               else {
                  (*phase)(i,j,e) = Util::MV;
               }
            }
         }
      }
   }
   return true;
}

float CalibratorPhase::getMinPrecip() const {
   return mMinPrecip;
}
void CalibratorPhase::setMinPrecip(float iMinPrecip) {
   mMinPrecip = iMinPrecip;
}
void CalibratorPhase::setUseWetbulb(bool iUseWetbulb) {
   mUseWetbulb = iUseWetbulb;
}
bool CalibratorPhase::getUseWetbulb() {
   return mUseWetbulb;
}

std::string CalibratorPhase::description() {
   std::stringstream ss;
   ss << "   -c phase                     Compute precipitation phase based on temperature, with" << std::endl;
   ss << "                                values encoded by:" << std::endl;
   ss << "                                * 0 = no precipitation (Precip < minPrecip)" << std::endl;
   ss << "                                * 1 = rain (b < T)" << std::endl;
   ss << "                                * 2 = sleet (a < T <= b)" << std::endl;
   ss << "                                * 3 = snow (T < a" << std::endl;
   ss << "                                T can be either regular temperature or wetbulb temperature." << std::endl;
   ss << "                                Precip, and Temperature must be available to determine phase. If" << std::endl;
   ss << "                                using wetbulb, then relative humidity must also be available." << std::endl;
   ss << "                                Pressure is currently not needed because a standard atmosphere is used." << std::endl;
   ss << "      parameters=required       Read parameters from this text file. The file format is:" << std::endl;
   ss << "                                offset0 a b" << std::endl;
   ss << "                                    ...    " << std::endl;
   ss << "                                offsetN a b" << std::endl;
   ss << "                                If the file only has a single line, then the same set of parameters" << std::endl;
   ss << "                                are used for all offsets.                                          " << std::endl;
   ss << "      useWetbulb=1              If 1 use the wetbulb temperature to determine phase. If 0 use regular" << std::endl;
   ss << "                                temperature." << std::endl;
   ss << "      minPrecip=0.2             Minimum precip (in mm) needed to be considered as precipitation." << std::endl;
   return ss.str();
}
float CalibratorPhase::getWetbulb(float iTemperature, float iPressure, float iRelativeHumidity) {
   float temperatureC = iTemperature - 273.15;
   if(temperatureC <= -243.04 || iRelativeHumidity <= 0)
      return Util::MV;
   if(Util::isValid(temperatureC) && Util::isValid(iPressure) && Util::isValid(iRelativeHumidity)) {
      float e  = (iRelativeHumidity)*0.611*exp((17.63*temperatureC)/(temperatureC+243.04));
      float Td = (116.9 + 243.04*log(e))/(16.78-log(e));
      float gamma = 0.00066 * iPressure/1000;
      float delta = (4098*e)/pow(Td+243.04,2);
      if(gamma + delta == 0)
         return Util::MV;
      float wetbulbTemperature   = (gamma * temperatureC + delta * Td)/(gamma + delta);
      float wetbulbTemperatureK  = wetbulbTemperature + 273.15;
      return wetbulbTemperatureK;
   }
   else {
      return Util::MV;
   }
}
