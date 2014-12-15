#include "Coeffs.h"
#include <iostream>

Coeffs::Coeffs(Site iSite, std::vector<float> iCoeffs) : mSite(iSite), mCoeffs(iCoeffs) {}

float Coeffs::computeFactor(float iDirection, float iSpeed) const {
   float factor = mCoeffs[0];
   float pi     = 3.14159265;
   float rads   = iDirection * pi / 180;
   int i = 1;
   float freq = 1;
   while(i < mCoeffs.size()) {
      factor += mCoeffs[i] * sin(rads * freq);
      i++;
      factor += mCoeffs[i] * cos(rads * freq);
      i++;
      freq   += 1;
   }
   if(factor > mMaxFactor) {
      std::cout << "Warning: A factor of " << factor << " for site " << mSite.getId()
                << " is too high. Resetting factor to ";
      factor = mMaxFactor;
      std::cout << factor << "." << std::endl;
   }
   if(factor < mMinFactor) {
      std::cout << "Warning: A factor of " << factor << " for site " << mSite.getId()
                << " is too small. Resetting factor to ";
      factor = mMinFactor;
      std::cout << factor << "." << std::endl;
   }
   float increment = iSpeed * (factor - 1);
   if(increment > mMaxSpeedIncrease) {
      std::cout << "Warning: A speed increase of " << increment << " for site " << mSite.getId()
                << " is > " << mMaxSpeedIncrease << ". Resetting factor to ";
      factor = 1 + mMaxSpeedIncrease / iSpeed;
      std::cout << factor << "." << std::endl;
   }
   if(-increment > mMaxSpeedDecrease) {
      std::cout << "Warning: A speed decrease of " << -increment << " for site " << mSite.getId()
                << " is > " << fabs(mMaxSpeedDecrease) << ". Resetting factor to ";
      factor = 1 - mMaxSpeedDecrease / iSpeed;
      std::cout << factor << "." << std::endl;
   }
   return factor;
}
