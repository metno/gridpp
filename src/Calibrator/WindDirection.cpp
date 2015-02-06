#include "WindDirection.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile.h"
#include "../Parameters.h"
CalibratorWindDirection::CalibratorWindDirection(const ParameterFile* iParameterFile, Variable::Type iVariable):
      Calibrator(),
      mVariable(iVariable),
      mParameterFile(iParameterFile) {

   if(mParameterFile->getNumParameters() != 9) {
      Util::error("CalibratorWindDirection: ParameterFile must have 9 parameters");
   }
}

bool CalibratorWindDirection::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& wind      = *iFile.getField(mVariable, t);
      Field& direction = *iFile.getField(Variable::WD, t);

      Parameters parameters = mParameterFile->getParameters(t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               float currDirection = direction(i,j,e);
               float factor = getFactor(currDirection, parameters);
               wind(i,j,e) = factor*wind(i,j,e);
            }
         }
      }
   }
   return true;
}

std::string CalibratorWindDirection::description() {
   std::stringstream ss;
   ss << "   -c windDirection             Multiply a variable by a factor based on the wind-direction:" << std::endl;
   ss << "                                factor = a + b*sin(dir)   + c*cos(dir)   + d*sin(2*dir) + e*cos(2*dir)" << std::endl;
   ss << "                                           + f*sin(3*dir) + g*cos(3*dir) + h*sin(4*dir) + i*cos(4*dir)" << std::endl;
   ss << "      parameters=required       Read parameters from this text file. The file format is:" << std::endl;
   ss << "                                offset0 a b c d e f g h i" << std::endl;
   ss << "                                           ...          " << std::endl;
   ss << "                                offsetN a b c d e f g h i" << std::endl;
   ss << "                                If the file only has a single line, then the same set of parameters" << std::endl;
   ss << "                                are used for all offsets.                                          " << std::endl;
   return ss.str();
}

float CalibratorWindDirection::getFactor(float iWindDirection, const Parameters& iPar) {
   if(!Util::isValid(iWindDirection) || !iPar.isValid())
      return Util::MV;
   float dir = Util::deg2rad(iWindDirection);
   float factor = 1;
   factor = iPar[0] +
            iPar[1] * sin(dir) +
            iPar[2] * cos(dir) +
            iPar[3] * sin(2*dir) +
            iPar[4] * cos(2*dir) +
            iPar[5] * sin(3*dir) +
            iPar[6] * cos(3*dir) +
            iPar[7] * sin(4*dir) +
            iPar[8] * cos(4*dir);
   if(factor < 0)
      factor = 0;
   return factor;
}
