#include "Phase.h"
#include "../Util.h"
#include "../File/File.h"
#include <cmath>
CalibratorPhase::CalibratorPhase(const ParameterFile* iParameterFile) :
      Calibrator(),
      mParameterFile(iParameterFile) {

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


   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      const Parameters& par = mParameterFile->getParameters(t);
      float snowSleetThreshold = par[0];
      float sleetRainThreshold = par[1];
      const Field& precip   = *iFile.getField(Variable::Precip, t);
      const Field& temp     = *iFile.getField(Variable::T, t);
      const Field& pressure = *iFile.getField(Variable::P, t);
      const Field& rh       = *iFile.getField(Variable::RH, t);
      Field& phase          = *iFile.getField(Variable::Phase, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               float currPrecip   = precip(i,j,e);
               float currTemp     = temp(i,j,e);
               float currPressure = pressure(i,j,e);
               float currRh       = rh(i,j,e);
               if(Util::isValid(snowSleetThreshold) && Util::isValid(sleetRainThreshold) && Util::isValid(currPrecip) && Util::isValid(currTemp) && Util::isValid(currPressure) && Util::isValid(currRh)) {
                  float currWetbulb = getWetbulb(currTemp, currPressure, currRh);
                  if(currPrecip == 0)
                     phase(i,j,e)  = CalibratorPhase::PhaseNone;
                  else if(!Util::isValid(currWetbulb))
                     phase(i,j,e)  = Util::MV;
                  else if(currWetbulb <= snowSleetThreshold)
                     phase(i,j,e)  = CalibratorPhase::PhaseSnow;
                  else if(currWetbulb <= sleetRainThreshold)
                     phase(i,j,e)  = CalibratorPhase::PhaseSleet;
                  else
                     phase(i,j,e)  = CalibratorPhase::PhaseRain;
               }
               else {
                  phase(i,j,e) = Util::MV;
               }
            }
         }
      }
   }
   return true;
}
std::string CalibratorPhase::description() {
   std::stringstream ss;
   ss << "   -c phase                     Compute precipitation phase based on wetbulb temperature, with" << std::endl;
   ss << "                                values encoded by:" << std::endl;
   ss << "                                * 0 = no precipitation (Precip == 0)" << std::endl;
   ss << "                                * 1 = rain (b < Tw)" << std::endl;
   ss << "                                * 2 = sleet (a < Tw <= b)" << std::endl;
   ss << "                                * 3 = snow (Tw < a" << std::endl;
   ss << "                                Precip amount, temperature, relative humidity, and pressure must be available" << std::endl;
   ss << "      parameters=required       Read parameters from this text file. The file format is:" << std::endl;
   ss << "                                offset0 a b" << std::endl;
   ss << "                                    ...    " << std::endl;
   ss << "                                offsetN a b" << std::endl;
   ss << "                                If the file only has a single line, then the same set of parameters" << std::endl;
   ss << "                                are used for all offsets.                                          " << std::endl;
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
