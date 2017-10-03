#ifndef DOWNSCALER_BILINEAR_H
#define DOWNSCALER_BILINEAR_H
#include <map>
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
class File;
class DownscalerBilinear : public Downscaler {
   public:
      DownscalerBilinear(Variable::Type iVariable, const Options& iOptions);
      static std::string description();
      std::string name() const {return "bilinear";};
      //! Bilinearly interpolate four irregularly placed points:
      //! v1 v3
      //! v0 v2
      //! with corresponding x,y coordinates.
      //! @param x Interpolate to this x-coordinate
      //! @param y Interpolate to this y-coordinate
      //! Two or more points cannot be colocated
      static float bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3);
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
};
#endif
