#include "Variable.h"
#include "Util.h"
#include <sstream>

std::string Variable::getTypeName(Type iType) {
   if(iType == Precip)
      return "Precip";
   else if(iType == PrecipAcc)
      return "PrecipAcc";
   else if(iType == Pop)
      return "Pop";
   else if(iType == Pop6h)
      return "Pop6h";
   else if(iType == PrecipLow)
      return "PrecipLow";
   else if(iType == PrecipMiddle)
      return "PrecipMiddle";
   else if(iType == PrecipHigh)
      return "PrecipHigh";
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
   else if(iType == MSLP)
      return "MSLP";
   else if(iType == QNH)
      return "QNH";
   else if(iType == LwinAcc)
      return "LwinAcc";
   else if(iType == SwinAcc)
      return "SwinAcc";
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
   else if(iName == "Pop")
      return Pop;
   else if(iName == "Pop6h")
      return Pop6h;
   else if(iName == "PrecipLow")
      return PrecipLow;
   else if(iName == "PrecipMiddle")
      return PrecipMiddle;
   else if(iName == "PrecipHigh")
      return PrecipHigh;
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
   else if(iName == "MSLP")
      return MSLP;
   else if(iName == "QNH")
      return QNH;
   else if(iName == "SwinAcc")
      return SwinAcc;
   else if(iName == "LwinAcc")
      return LwinAcc;
   else if(iName == "Fake")
      return Fake;
   else
      Util::error("Undefined variable type: " + iName);
}

std::string Variable::getDescriptions() {
   std::stringstream ss;
   ss << Util::formatDescription("-v T", "Temperature") << std::endl;
   ss << Util::formatDescription("-v Precip", "Hourly precip") << std::endl;
   ss << Util::formatDescription("-v PrecipAcc", "Accumulated precip") << std::endl;
   ss << Util::formatDescription("-v Pop", "Probability of precip") << std::endl;
   ss << Util::formatDescription("-v Pop6h", "Probability of precip for the last 6 hours") << std::endl;
   ss << Util::formatDescription("-v PrecipLow", "Low estimate of precipitation") << std::endl;
   ss << Util::formatDescription("-v PrecipMiddl", "Middle estimate of precipitation") << std::endl;
   ss << Util::formatDescription("-v PrecipHigh", "High estimate of precipitation") << std::endl;
   ss << Util::formatDescription("-v W", "Wind speed") << std::endl;
   ss << Util::formatDescription("-v WD", "Wind direction") << std::endl;
   ss << Util::formatDescription("-v U", "U-wind") << std::endl;
   ss << Util::formatDescription("-v V", "V-wind") << std::endl;
   ss << Util::formatDescription("-v Cloud", "Cloud cover") << std::endl;
   ss << Util::formatDescription("-v RH", "Relative humidity") << std::endl;
   ss << Util::formatDescription("-v Phase", "Precipitation phase (0 none, 1 rain, 2 sleet, 3 snow)") << std::endl;
   ss << Util::formatDescription("-v P", "Pressure") << std::endl;
   ss << Util::formatDescription("-v MSLP", "Mean sea-level pressure") << std::endl;
   ss << Util::formatDescription("-v QNH", "Pressure reduced to sea-level using standard atmosphere (ICAO)") << std::endl;
   ss << Util::formatDescription("-v SwinAcc", "Accumulated incoming shortwave radiation") << std::endl;
   ss << Util::formatDescription("-v LwinAcc", "Accumulated incoming longwave radiation") << std::endl;
   return ss.str();
}

std::vector<Variable::Type> Variable::getAllVariables() {
   std::vector<Type> variables;
   variables.push_back(Variable::T);
   variables.push_back(Variable::PrecipAcc);
   variables.push_back(Variable::Pop);
   variables.push_back(Variable::Pop6h);
   variables.push_back(Variable::PrecipLow);
   variables.push_back(Variable::PrecipMiddle);
   variables.push_back(Variable::PrecipHigh);
   variables.push_back(Variable::Precip);
   variables.push_back(Variable::SwinAcc);
   variables.push_back(Variable::LwinAcc);
   variables.push_back(Variable::W);
   variables.push_back(Variable::WD);
   variables.push_back(Variable::U);
   variables.push_back(Variable::V);
   variables.push_back(Variable::Cloud);
   variables.push_back(Variable::RH);
   variables.push_back(Variable::Phase);
   variables.push_back(Variable::P);
   variables.push_back(Variable::MSLP);
   variables.push_back(Variable::QNH);
   variables.push_back(Variable::LwinAcc);
   variables.push_back(Variable::SwinAcc);
   // Do not include Variable::NONE;
   return variables;
}

float Variable::getMin(Type iType) {
   switch(iType) {
      case T:
         return 0;
      case Precip:
         return 0;
      case PrecipAcc:
         return 0;
      case Pop:
         return 0;
      case Pop6h:
         return 0;
      case PrecipLow:
         return 0;
      case PrecipMiddle:
         return 0;
      case PrecipHigh:
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
      case MSLP:
         return 0;
      case QNH:
         return 0;
      case LwinAcc:
         return Util::MV;
      case SwinAcc:
         return Util::MV;
      default:
         return Util::MV;
   }
}
float Variable::getMax(Type iType) {

   switch(iType) {
      case T:
         return Util::MV;
      case Precip:
         return Util::MV;
      case PrecipAcc:
         return Util::MV;
      case Pop:
         return 1;
      case Pop6h:
         return 1;
      case PrecipLow:
         return Util::MV;
      case PrecipMiddle:
         return Util::MV;
      case PrecipHigh:
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
      case MSLP:
         return Util::MV;
      case QNH:
         return Util::MV;
      case LwinAcc:
         return Util::MV;
      case SwinAcc:
         return Util::MV;
      default:
         return Util::MV;
   }
}
std::string Variable::getUnits(Type iType) {
   switch(iType) {
      case T:
         return "K";
      case Precip:
         return "kg/m^2";
      case PrecipAcc:
         return "kg/m^2";
      case Pop:
         return "1";
      case Pop6h:
         return "1";
      case PrecipLow:
         return "kg/m^2";
      case PrecipMiddle:
         return "kg/m^2";
      case PrecipHigh:
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
         return "1";
      case Phase:
         return "";
      case P:
         return "pa";
      case MSLP:
         return "pa";
      case QNH:
         return "pa";
      case LwinAcc:
         return "W s/m^2";
      case SwinAcc:
         return "W s/m^2";
      default:
         return "no_units";
   }
}
std::string Variable::getStandardName(Type iType) {
   switch(iType) {
      case T:
         return "air_temperature";
      case Precip:
         return "precipitation_amount";
      case PrecipAcc:
         return "precipitation_amount_acc";
      case Pop:
         return "precipitation_amount";
      case Pop6h:
         return "precipitation_amount";
      case PrecipLow:
         return "precipitation_amount";
      case PrecipMiddle:
         return "precipitation_amount";
      case PrecipHigh:
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
         return "";
      case P:
         return "surface_air_pressure";
      case MSLP:
         return "air_pressure_at_sea_level";
      case QNH:
         // TODO: What is the standard name for QNH?
         return "air_pressure_at_sea_level_qnh";
      case SwinAcc:
         return "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time";
      case LwinAcc:
         return "integral_of_surface_downwelling_longwave_flux_in_air_wrt_time";
      default:
         return "no_units";
   }
}
