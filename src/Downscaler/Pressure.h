#ifndef DOWNSCALER_PRESSURE_H
#define DOWNSCALER_PRESSURE_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
//! Adjust the value of the nearest neighbour (nn), based on a standard atmosphere and the elevation
//! difference to the lookup point (p):
//! P(p) = P(nn) * exp(-1.21e-4*(elev(p) - elev(nn))
//! Uses the pressure at the nearest neighbour when the lookup location does not have an elevation
class DownscalerPressure : public Downscaler {
   public:
      DownscalerPressure(Variable::Type iVariable);
      static std::string description();
      std::string name() const {return "pressure";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      static const float mConstant;
};
#endif
