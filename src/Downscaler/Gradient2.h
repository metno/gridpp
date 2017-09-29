#ifndef DOWNSCALER_GRADIENT2_H
#define DOWNSCALER_GRADIENT2_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
#include "../Field.h"
typedef std::vector<std::vector<int> > vec2Int;
//! Adjust the value of the nearest neighbour (nn), based on the gradient in a neighbourhood
//! surrounding the nearest neighbour, and the elevation difference to the lookup point (p):
//! T(p) = T(nn) + gradient*(elev(p) - elev(nn))
//! If the variable is log transformed, then use:
//! T(p) = T(nn) * exp(gradient*(elev(p) - elev(nn)))
//! Uses nearest neighbour when the lookup location does not have an elevation
class DownscalerGradient2 : public Downscaler {
   public:
      //! Downscale the specified variable
      DownscalerGradient2(Variable::Type iVariable, const Options& iOptions);
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
      float mMinLafForElevGradient;
      int mMinNumPoints;
      float mMinFracSeaPoints;
      float mDefaultGradient;
      int mLafRadius;
      std::vector<int>  mLafSearchRadii;
      std::vector<float>  mLafWeights;
      float mMinLafDiff; // Minimum elevation difference within neighbourhood to use gradient
      bool mAverageNeighbourhood;
      mutable bool mHasIssuedWarningUnstable;
      bool mSaveGradient;
      Variable::Type mGradientVariable;
      float calcLafGradient(int i, int j, int e, int Icenter, int Jcenter, const Field& iField, const vec2& iLafs, const vec2& iElevs, float iElevGradient) const;
      float calcElevGradient(int i, int j, int e, int Icenter, int Jcenter, const Field& iField, const Field& iGfield, const vec2& iElevs, const vec2& iLafs) const;
      bool calcNeighbourhoodMean(const Field& iField, const vec2& iElevs, int i, int j, int e, int Icenter, int Jcenter, int iRadius, float& iValueMean, float& iElevMean) const;
      bool calcBaseValues(const Field& iField, const vec2& iElevs, const vec2& iLafs, int i, int j, int e, int Icenter, int Jcenter, int iRadius, float& iValueMean, float& iElevMean, float& iLafMean) const;
      // float calcLafGradient(float iLaf, int iIcenter, int iJcenter, int e, vec2 iLafs, vec2 iElevs, const Field& iField, const File& iInput) const;
      // float calcElevGradient(float iLaf, int iIcenter, int iJcenter, int e, vec2 iLafs, vec2 iElevs, const Field& iField, const File& iInput) const;
};
#endif
