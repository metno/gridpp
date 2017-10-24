#ifndef DOWNSCALER_BILINEAR_H
#define DOWNSCALER_BILINEAR_H
#include <map>
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
#include "../Field.h"
class File;
class DownscalerBilinear : public Downscaler {
   public:
      DownscalerBilinear(Variable::Type iVariable, const Options& iOptions);

      //! Bilinearly interpolate four irregularly placed points:
      //! v1 v3
      //! v0 v2
      //! with corresponding x,y coordinates.
      //! @param x Interpolate to this x-coordinate
      //! @param y Interpolate to this y-coordinate
      //! Two or more points cannot be colocated
      static float bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3);

      // Interpolate a whole field
      static void downscaleField(const Field& iInput, Field& iOutput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ);

      // Interpolate a whole vec2
      static vec2 downscaleVec(const vec2& iInput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ);

      static std::string description();
      std::string name() const {return "bilinear";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
};
#endif
