#ifndef DOWNSCALER_GRADIENT_H
#define DOWNSCALER_GRADIENT_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
//! Adjust the value of the nearest neighbour (nn), based on the gradient in a neighbourhood
//! surrounding the nearest neighbour, and the elevation difference to the lookup point (p):
//! T(p) = T(nn) + gradient*(elev(p) - elev(nn))
class DownscalerGradient : public Downscaler {
   public:
      //! Downscale the specified variable
      DownscalerGradient(Variable::Type iVariable);
      //! Do not compute the gradient but set it to a fixed amount. Positive rate means
      //! increasing with elevation. The units are per meters.
      void setConstantGradient(float iGradient);
      float getConstantGradient() const;
      //! Calculate gradient in a neighbourhood of points within +- iNumPoints
      //! in each direction.
      void setSearchRadius(int iNumPoints);
      int  getSearchRadius() const;
      static std::string description();
      std::string name() const {return "gradient";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      int mSearchRadius;
      float mConstGradient;
};
#endif
