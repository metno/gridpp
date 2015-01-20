#include "Variable.h"
#include "Util.h"

std::string Variable::getTypeName(Type iType) {
   if(iType == Precip)
      return "Precip";
   else if(iType == PrecipAcc)
      return "Accumulated precip";
   else if(iType == Cloud)
      return "Cloud cover";
   else if(iType == T)
      return "Temperature";
   else if(iType == U)
      return "U-wind";
   else if(iType == V)
      return "V-wind";
   else if(iType == W)
      return "Wind speed";
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
      Util::error("Undefined variable type");
}
