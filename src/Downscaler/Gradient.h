#ifndef DOWNSCALER_GRADIENT_H
#define DOWNSCALER_GRADIENT_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
//! Adjust the value of the nearest neighbour (nn), based on the gradient in a neighbourhood
//! surrounding the nearest neighbour, and the elevation difference to the lookup point (p):
//! T(p) = T(nn) + gradient*(elev(p) - elev(nn))
//! If the variable is log transformed, then use:
//! T(p) = T(nn) * exp(gradient*(elev(p) - elev(nn)))
//! Uses nearest neighbour when the lookup location does not have an elevation
class DownscalerGradient : public Downscaler {
   public:
      //! Downscale the specified variable
      DownscalerGradient(Variable::Type iVariable, const Options& iOptions);
      float getConstantGradient() const;
      int   getSearchRadius() const;
      float getMinElevDiff() const;
      bool  getLogTransform() const;
      float getMinGradient() const;
      float getMaxGradient() const;
      float getDefaultGradient() const;
      static std::string description();
      std::string name() const {return "gradient";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      int   mSearchRadius;
      float mConstantGradient;
      float mMinElevDiff; // Minimum elevation difference within neighbourhood to use gradient
      bool  mLogTransform;
      float mMinGradient;
      float mMaxGradient;
      float mDefaultGradient;
      bool mAverageNeighbourhood;
      mutable bool mHasIssuedWarningUnstable;
      std::string mDownscalerName;
};
#endif
