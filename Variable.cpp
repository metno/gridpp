#include "Variable.h"

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
