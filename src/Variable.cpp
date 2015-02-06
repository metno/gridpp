#include "Variable.h"
#include "Util.h"
#include <sstream>

std::string Variable::getTypeName(Type iType) {
   if(iType == Precip)
      return "Precip";
   else if(iType == PrecipAcc)
      return "PrecipAcc";
   else if(iType == Cloud)
      return "Cloud";
   else if(iType == T)
      return "T";
   else if(iType == U)
      return "U";
   else if(iType == V)
      return "V";
   else if(iType == W)
      return "W";
   else if(iType == WD)
      return "WD";
   else if(iType == RH)
      return "RH";
   else if(iType == Phase)
      return "Phase";
   else if(iType == P)
      return "P";
   else if(iType == Fake)
      return "Fake";
   else
      return "Unknown";
}

Variable::Type Variable::getType(std::string iName) {
   if(iName == "Precip")
      return Precip;
   else if(iName == "PrecipAcc")
      return PrecipAcc;
   else if(iName == "Cloud")
      return Cloud;
   else if(iName == "T")
      return T;
   else if(iName == "U")
      return U;
   else if(iName == "V")
      return V;
   else if(iName == "W")
      return W;
   else if(iName == "WD")
      return WD;
   else if(iName == "RH")
      return RH;
   else if(iName == "Phase")
      return Phase;
   else if(iName == "P")
      return P;
   else if(iName == "Fake")
      return Fake;
   else
      Util::error("Undefined variable type: " + iName);
}

std::string Variable::description() {
   std::stringstream ss;
   ss << "   -v T                         Temperature" << std::endl;
   ss << "   -v Precip                    Hourly precip" << std::endl;
   ss << "   -v PrecipAcc                 Accumulated precip" << std::endl;
   ss << "   -v W                         Wind speed" << std::endl;
   ss << "   -v WD                        Wind direction" << std::endl;
   ss << "   -v U                         U-wind" << std::endl;
   ss << "   -v V                         V-wind" << std::endl;
   ss << "   -v Cloud                     Cloud cover" << std::endl;
   ss << "   -v RH                        Relative humidity" << std::endl;
   ss << "   -v Phase                     Precipitation phase (0 none, 1 rain, 2 sleet, 3 snow)" << std::endl;
   ss << "   -v P                         Pressure" << std::endl;
   return ss.str();
}

std::vector<Variable::Type> Variable::getAllVariables() {
   std::vector<Type> variables;
   variables.push_back(Variable::T);
   variables.push_back(Variable::PrecipAcc);
   variables.push_back(Variable::Precip);
   variables.push_back(Variable::W);
   variables.push_back(Variable::WD);
   variables.push_back(Variable::U);
   variables.push_back(Variable::V);
   variables.push_back(Variable::Cloud);
   variables.push_back(Variable::RH);
   variables.push_back(Variable::Phase);
   variables.push_back(Variable::P);
   return variables;
}

float Variable::getMin(Type iType) {
   switch(iType) {
      case T:
         return 0;
      case PrecipAcc:
         return 0;
      case Precip:
         return 0;
      case W:
         return 0;
      case WD:
         return 0;
      case U:
         return Util::MV;
      case V:
         return Util::MV;
      case Cloud:
         return 0;
      case RH:
         return 0;
      case Phase:
         return 0;
      case P:
         return 0;
      default:
         return Util::MV;
   }
}
float Variable::getMax(Type iType) {

   switch(iType) {
      case T:
         return Util::MV;
      case PrecipAcc:
         return Util::MV;
      case Precip:
         return Util::MV;
      case W:
         return Util::MV;
      case WD:
         return 360;
      case U:
         return Util::MV;
      case V:
         return Util::MV;
      case Cloud:
         return 1;
      case RH:
         return 1.1;
      case Phase:
         return 3;
      case P:
         return Util::MV;
      default:
         return Util::MV;
   }
}
std::string Variable::getUnits(Type iType) {
   switch(iType) {
      case T:
         return "K";
      case PrecipAcc:
         return "kg/m^2";
      case Precip:
         return "kg/m^2";
      case W:
         return "m/s";
      case WD:
         return "degrees";
      case U:
         return "m/s";
      case V:
         return "m/s";
      case Cloud:
         return "1";
      case RH:
         return "%";
      case Phase:
         return "phase";
      case P:
         return "pa";
      default:
         return "no_units";
   }
}
std::string Variable::getStandardName(Type iType) {
   switch(iType) {
      case T:
         return "air_temperature";
      case PrecipAcc:
         return "precipitation_amount_acc";
      case Precip:
         return "precipitation_amount";
      case W:
         return "wind_speed";
      case WD:
         return "wind_from_direction";
      case U:
         // TODO: Is this correct?
         return "eastward_wind";
      case V:
         // TODO: Is this correct?
         return "northward_wind";
      case Cloud:
         return "cloud_area_fraction";
      case RH:
         return "relative_humidity";
      case Phase:
         return "phase";
      case P:
         return "surface_air_pressure";
      default:
         return "no_units";
   }
}
