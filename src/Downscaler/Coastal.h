#ifndef DOWNSCALER_COASTAL_H
#define DOWNSCALER_COASTAL_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
class DownscalerCoastal : public Downscaler {
   public:
      //! Downscale the specified variable
      DownscalerCoastal(Variable::Type iVariable, const Options& iOptions);
      int   getSearchRadius() const;
      float getMinLafDiff() const;
      float getMinGradient() const;
      float getMaxGradient() const;
      static std::string description();
      std::string name() const {return "coastal";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      int   mSearchRadius;
      float mMinLafDiff; // Minimum elevation difference within neighbourhood to use gradient
      float mMinGradient;
      float mMaxGradient;
};
#endif
