#include "Pressure.h"
#include "../File/File.h"
#include "../Util.h"
#include "gridpp.h"
#include <math.h>

DownscalerPressure::DownscalerPressure(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions),
      mTemperatureVariable("") {
   iOptions.getValue("temperatureVariable", mTemperatureVariable);
   iOptions.check();
}

void DownscalerPressure::downscaleCore(const File& iInput, File& iOutput) const {
   int nY = iOutput.getNumY();
   int nX = iOutput.getNumX();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();
   vec2 ielevs = iInput.getElevs();
   vec2 oelevs = iOutput.getElevs();

   gridpp::Grid igrid(iInput.getLats(), iInput.getLons());
   gridpp::Grid ogrid(iOutput.getLats(), iOutput.getLons());

   // Interpolate elevations to new grid
   vec2 oelevsInterp = gridpp::nearest(igrid, ogrid, ielevs);
   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);
      for(int e = 0; e < nEns; e++) {
         vec2 ofieldInterp = gridpp::nearest(igrid, ogrid, ifield(e));
         vec2 tfieldInterp;
         if(mTemperatureVariable != "") {
             Field& tfield = *iOutput.getField(mTemperatureVariable, t);
             vec2 tfieldInterp = gridpp::nearest(igrid, ogrid, tfield(e));
         }
         for(int y = 0; y < nY; y++) {
            for(int x = 0; x < nX; x++) {
                float temperature = 288.15;
                if(tfieldInterp.size() != 0)
                    temperature = tfieldInterp[y][x];
                ofield(y, x, e) = gridpp::pressure(oelevsInterp[y][x], oelevs[y][x], ofieldInterp[y][x], temperature);
            }
         }
      }
   }
}
std::string DownscalerPressure::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-d pressure", "Adjusts the pressure of the nearest neighbour based on the elevation difference and a standard atmosphere.") << std::endl;
      ss << Util::formatDescription("   temperatureVariable=undef", "Which variable to read temperature from? If undefined, a temperature of 288.15 K is used.") << std::endl;
   }
   else
      ss << Util::formatDescription("-d pressure", "Adjusts the pressure based on the elevation differences") << std::endl;
   return ss.str();
}
