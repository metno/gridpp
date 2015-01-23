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
   return variables;
}
