#ifndef COEFFS_H
#define COEFFS_H
#include <vector>
#include <math.h>
#include "Site.h"

//! Represents Fourier coefficients for a site
class Coeffs {
   public:
      Coeffs(Site iSite, std::vector<float> iCoeffs);
      //! Compute the multiplication factor for a given direction and wind speed
      //! iDirection: in degrees
      //! Note: Speed is only used to check that the speed increment is within bounds.
      float computeFactor(float iDirection, float iSpeed) const;
      Site getSite() const {return mSite;};
      //! Don't allow a factor greater than this:
      float getMaxFactor() const {return mMaxFactor;};
      //! Don't allow a factor smaller than this:
      float getMinFactor() const {return mMinFactor;};
      //! Don't allow a speed increase (m/s) greater than this:
      float getMaxSpeedIncrease() const {return mMaxSpeedIncrease;};
      //! Don't allow a speed decrease (m/s) greater than this:
      float getMaxSpeedDecrease() const {return mMaxSpeedDecrease;};
   private:
      Site mSite;
      std::vector<float> mCoeffs;
      static const float mMaxFactor = 3;
      static const float mMinFactor = 0.25;
      static const float mMaxSpeedIncrease = 8;
      static const float mMaxSpeedDecrease = 8;
};
#endif
