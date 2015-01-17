#ifndef DOWNSCALER_GRADIENT_H
#define DOWNSCALER_GRADIENT_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
class DownscalerGradient : public Downscaler {
   public:
      DownscalerGradient(Variable::Type iVariable);
      void setConstantGradient(float iGradient);
      void setNeighbourhoodSize(int iSize);
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      Variable::Type mVariable;
      int mSize;
      float mConstGradient;
};
#endif
