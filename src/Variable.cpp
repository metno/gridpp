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
   ss << "   -v W                         Wind speed" << std::endl;
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
   variables.push_back(Variable::W);
   variables.push_back(Variable::U);
   variables.push_back(Variable::V);
   variables.push_back(Variable::Cloud);
   variables.push_back(Variable::RH);
   variables.push_back(Variable::Phase);
   variables.push_back(Variable::P);
   return variables;
}
