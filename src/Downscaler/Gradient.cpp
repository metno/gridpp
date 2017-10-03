#include "Gradient.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerGradient2::DownscalerGradient2(Variable::Type iVariable, const Options& iOptions) :
      Downscaler(iVariable, iOptions),
      mSearchRadius(3),
      mMinGradient(Util::MV),
      mMaxGradient(Util::MV),
      mLogTransform(false),
      mConstantGradient(Util::MV),
      mDefaultGradient(0),
      mMinElevDiff(30),
      mAverageNeighbourhood(false),
      mSaveGradient(""),
      mMinLafForElevGradient(0),
      mMinLafDiff(0.1),
      mLafRadius(1),
      mMinNumPoints(2),
      mMinFracSeaPoints(0),
      mGradientVariable(Variable::None),
      mHasIssuedWarningUnstable(false) {

   iOptions.getValue("minGradient", mMinGradient);
   iOptions.getValue("maxGradient", mMaxGradient);
   iOptions.getValue("defaultGradient", mDefaultGradient);
   iOptions.getValue("searchRadius", mSearchRadius);
   iOptions.getValue("logTransform", mLogTransform);
   iOptions.getValue("minLafDiff", mMinLafDiff);
   iOptions.getValue("lafRadius", mLafRadius);
   iOptions.getValue("minLafForElevGradient", mMinLafForElevGradient);
   iOptions.getValue("constantGradient", mConstantGradient);
   iOptions.getValue("minElevDiff", mMinElevDiff);
   iOptions.getValue("averageNeighbourhood", mAverageNeighbourhood);
   iOptions.getValue("saveGradient", mSaveGradient);
   iOptions.getValue("minNumPoints", mMinNumPoints);
   iOptions.getValue("minFracSeaPoints", mMinFracSeaPoints);
   std::string gradientVariable;
   if(iOptions.getValue("gradientVariable", gradientVariable)) {
      mGradientVariable = Variable::getType(gradientVariable);
   }
   else {
      mGradientVariable = mVariable;
   }

   iOptions.getValues("lafSearchRadii", mLafSearchRadii);
   iOptions.getValues("lafWeights", mLafWeights);
   // if(!iOptions.getValues("lafSearchRadii", mLafSearchRadii)) {
   //    mLafSearchRadii.push_back(1);
   //    mLafSearchRadii.push_back(2);
   //    mLafSearchRadii.push_back(3);
   // }

   // if(!iOptions.getValues("lafWeights", mLafWeights)) {
   //    mLafWeights.push_back(0.5);
   //    mLafWeights.push_back(0.3);
   //    mLafWeights.push_back(0.2);
   // }
   assert(mLafSearchRadii.size() == mLafWeights.size());
}

void DownscalerGradient2::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumLat();
   int nLon = iOutput.getNumLon();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   vec2 ilats  = iInput.getLats();
   vec2 ilons  = iInput.getLons();
   vec2 ielevs = iInput.getElevs();
   vec2 ilafs = iInput.getLandFractions();
   vec2 olats  = iOutput.getLats();
   vec2 olons  = iOutput.getLons();
   vec2 oelevs = iOutput.getElevs();
   vec2 olafs = iOutput.getLandFractions();

   float minAllowed = Variable::getMin(mVariable);
   float maxAllowed = Variable::getMax(mVariable);

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   getNearestNeighbour(iInput, iOutput, nearestI, nearestJ);
   if(!iInput.hasVariable(mGradientVariable)) {
      std::stringstream ss;
     ss <<"Cannot compute gradient, since file does not have " << Variable::getTypeName(mGradientVariable);
      Util::error(ss.str());
   }

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mVariable, t);
      Field& ofield = *iOutput.getField(mVariable, t);
      Field& gfield = *iInput.getField(mGradientVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int Icenter = nearestI[i][j];
            int Jcenter = nearestJ[i][j];
            assert(Icenter < ielevs.size());
            assert(Jcenter < ielevs[Icenter].size());
            for(int e = 0; e < nEns; e++) {
               float currElev = oelevs[i][j];
               float currLaf = olafs[i][j];
               float nearestElev = ielevs[Icenter][Jcenter];
               float nearestLaf = ilafs[Icenter][Jcenter];
               float nearestValue = ifield(Icenter, Jcenter, e);

               float elevGradient = calcElevGradient(i, j, e, Icenter, Jcenter, ifield, gfield, ielevs, ilafs);
               float lafGradient = Util::MV;
               if(mLafSearchRadii.size() > 0)
                  lafGradient = calcLafGradient(i, j, e, Icenter, Jcenter, ifield, ilafs, ielevs, elevGradient);
               if(mSaveGradient == "elev") {
                  ofield(i, j, e) = elevGradient;
               }
               else if(mSaveGradient == "laf") {
                  ofield(i, j, e) = lafGradient;
               }
               else {
                  /* General approach
                     1) Compute the elevation gradient, by ignoring any ocean points
                     2) Compute the LAF gradient, by adjusting neighbourhood using elevation gradient
                     3) Adjust the point with the lower LAF by adding on the LAF gradient and elevation gradient
                  */
                  // Compute elevation and value to extrapolate from
                  int averagingRadius = mAverageNeighbourhood * mSearchRadius;
                  float baseValue = Util::MV;
                  float baseElev = Util::MV;
                  float baseLaf = Util::MV;
                  // calcNeighbourhoodMean(ifield, ielevs, i, j, e, Icenter, Jcenter, averagingRadius, baseValue, baseElev);
                  baseValue = nearestValue;
                  baseElev = nearestElev;
                  baseLaf = nearestLaf;
                  float dElev = currElev - baseElev;
                  float dLaf = currLaf - baseLaf;

                  float value = baseValue;
                  // What base value do we use? For elevation gradient, we want to use the nearest
                  // neighbour, but for LAF gradient, we want to use the lowest LAF point.
                  value = baseValue + dElev * elevGradient;
                  if(Util::isValid(lafGradient))
                     value += dLaf * lafGradient;

                  if((Util::isValid(minAllowed) && value < minAllowed) || (Util::isValid(maxAllowed) && value > maxAllowed)) {
                     value = baseValue;
                  }
                  ofield(i,j,e)  = value;
               }
            }
         }
      }
   }
}

bool DownscalerGradient2::calcBaseValues(const Field& iField, const vec2& iElevs, const vec2& iLafs, int i, int j, int e, int Icenter, int Jcenter, int iRadius, float& iValueMean, float& iElevMean, float& iLafMean) const {
   float totalValue = 0;
   float totalElev = 0;
   float totalLaf = 0;
   int counter = 0;
   for(int ii = std::max(0, Icenter-iRadius); ii <= std::min(iField.getNumLat()-1, Icenter+iRadius); ii++) {
      for(int jj = std::max(0, Jcenter-iRadius); jj <= std::min(iField.getNumLon()-1, Jcenter+iRadius); jj++) {
         float currValue = iField(ii,jj,e);
         float currElev  = iElevs[ii][jj];
         if(Util::isValid(currValue) && Util::isValid(currElev)) {
            totalValue += currValue;
            totalElev  += currElev;
            counter++;
         }
      }
   }
   if(counter == 0) {
      iValueMean = Util::MV;
      iElevMean = Util::MV;
   }
   else {
      iValueMean = totalValue / counter;
      iElevMean = totalElev / counter;
   }

   return counter > 0;
}
bool DownscalerGradient2::calcNeighbourhoodMean(const Field& iField, const vec2& iElevs, int i, int j, int e, int Icenter, int Jcenter, int iRadius, float& iValueMean, float& iElevMean) const {
   float totalValue = 0;
   float totalElev = 0;
   int counter = 0;
   for(int ii = std::max(0, Icenter-iRadius); ii <= std::min(iField.getNumLat()-1, Icenter+iRadius); ii++) {
      for(int jj = std::max(0, Jcenter-iRadius); jj <= std::min(iField.getNumLon()-1, Jcenter+iRadius); jj++) {
         float currValue = iField(ii,jj,e);
         float currElev  = iElevs[ii][jj];
         if(Util::isValid(currValue) && Util::isValid(currElev)) {
            totalValue += currValue;
            totalElev  += currElev;
            counter++;
         }
      }
   }
   if(counter == 0) {
      iValueMean = Util::MV;
      iElevMean = Util::MV;
   }
   else {
      iValueMean = totalValue / counter;
      iElevMean = totalElev / counter;
   }

   return counter > 0;
}
float DownscalerGradient2::calcLafGradient(int i, int j, int e, int Icenter, int Jcenter, const Field& iField, const vec2& iLafs, const vec2& iElevs, float iElevGradient) const {
   /* Here is the formula:
      value = minT + (maxT  - minT) / (maxLaf - minLaf) * (laf - minLaf)
      value = minT + gradient * dLaf
   */
   // Compute the average temperature and LAF in a small radius
   float minAllowed = Variable::getMin(mVariable);
   float maxAllowed = Variable::getMax(mVariable);
   float minLaf = Util::MV;
   float maxLaf = Util::MV;
   float minElev = Util::MV;
   float maxElev = Util::MV;
   float minT = Util::MV;
   float maxT = Util::MV;
   // std::vector<float> values(mSearchRadii.size(), Util::MV);
   float value = 0;
   int counter = 0;
   float totalWeight = 0;
   for(int r = 0; r < mLafSearchRadii.size(); r++) {
      int numSeaPoints = 0;
      int numPoints = 0;
      int searchRadius = mLafSearchRadii[r];
      for(int ii = std::max(0, Icenter-searchRadius); ii <= std::min(iField.getNumLat()-1, Icenter+searchRadius); ii++) {
         for(int jj = std::max(0, Jcenter-searchRadius); jj <= std::min(iField.getNumLon()-1, Jcenter+searchRadius); jj++) {
            float currLaf  = iLafs[ii][jj];
            float currElev = iElevs[ii][jj];
            float currValue = iField(ii,jj,e);
            if(Util::isValid(currValue)) {
               numSeaPoints += currLaf < 0.5;
               numPoints++;
            }
            if(Util::isValid(currValue) && Util::isValid(currLaf)) {
               if(!Util::isValid(minLaf) || currLaf < minLaf) {
                  minLaf = currLaf;
                  minElev = currElev;
                  minT = currValue;
               }
               if(!Util::isValid(maxLaf) || currLaf > maxLaf) {
                  maxLaf = currLaf;
                  maxElev = currElev;
                  maxT = currValue;
               }
            }
         }
      }
      // std::cout << minLaf << " " << maxLaf << " " << minT << " " << maxT << std::endl;
      if(Util::isValid(maxLaf) && Util::isValid(minLaf) && (maxLaf - minLaf) >= mMinLafDiff && (float) numSeaPoints / numPoints >= mMinFracSeaPoints) {
         float lafGradient = (maxT - minT) / (maxLaf - minLaf);
         // Altitude-adjust maxT so that it is on the same elevation as minT
         // Before computing the gradient
         if(Util::isValid(maxElev) && Util::isValid(minElev)) {
            float maxTcorr = maxT + iElevGradient * (maxElev - minElev);
            lafGradient = (maxTcorr - minT) / (maxLaf - minLaf);
            value += lafGradient;
            totalWeight += mLafWeights[r];
            counter++;
         }
         /*
         // Safety check
         if(!Util::isValid(lafGradient))
            lafGradient = 0;
         // Check against minimum and maximum gradients
         if(Util::isValid(mMinGradient) && lafGradient < mMinGradient)
            lafGradient = mMinGradient;
         if(Util::isValid(mMaxGradient) && lafGradient > mMaxGradient)
            lafGradient = mMaxGradient;
         // std::cout << lafGradient << std::endl;
         // float currValue = nearestValue + dLaf * lafGradient;
         float currValue = minT + dLaf * lafGradient;

         // Altitude-adjust the current point
         if(Util::isValid(currElev) && Util::isValid(minElev)) {
            float dElev = currElev - minElev;
            currValue = currValue + dElev * mElevGradient;
         }
         // std::cout << currValue << std::endl;
         value += currValue * mWeights[r];
         totalWeight += mWeights[r];
         counter++;
         if(j == 43 && i == 56 && t == 0) {
            float dElev = currElev - minElev;
            float elevCorrection = dElev * mElevGradient;
            std::cout << minLaf << " " << maxLaf << " " << currLaf << "\n"
                      << minT << " " << maxT << " " << nearestValue << "\n"
                      << minElev << " " << maxElev << " " << currElev << "\n"
                      << currValue << " " << nearestValue << " " << value << "\n"
                      << lafGradient << " " << elevCorrection << std::endl;
         }
         */
      }
   }
   // std::cout << value << std::endl;
   float gradient;
   if(counter == 0 || !Util::isValid(value)) {
      gradient = Util::MV;
   }
   else {
      gradient = value / totalWeight;
   }
   return gradient;
}
float DownscalerGradient2::calcElevGradient(int i, int j, int e, int Icenter, int Jcenter, const Field& iField, const Field& iGfield, const vec2& iElevs, const vec2& iLafs) const {
   float elevGradient = Util::MV;
   if(Util::isValid(mConstantGradient)) {
      elevGradient = mConstantGradient;
   }
   else {
      /* Compute the model's gradient:
         The gradient is computed by using linear regression on forecast ~ elevation
         using all forecasts within a neighbourhood. To produce stable results, there
         is a requirement that the elevation within the neighbourhood has a large
         range (see mMinElevDiff).

         For bounded variables (e.g. wind speed), the gradient approach could cause
         forecasts to go outside its domain (e.g. negative winds). If this occurs,
         the nearest neighbour is used.
      */
      float meanXY  = 0; // elev*T
      float meanX   = 0; // elev
      float meanY   = 0; // T
      float meanXX  = 0; // elev*elev
      int   counter = 0;
      float min = Util::MV;
      float max = Util::MV;
      for(int ii = std::max(0, Icenter-mSearchRadius); ii <= std::min(iField.getNumLat()-1, Icenter+mSearchRadius); ii++) {
         for(int jj = std::max(0, Jcenter-mSearchRadius); jj <= std::min(iField.getNumLon()-1, Jcenter+mSearchRadius); jj++) {
            assert(ii < iElevs.size());
            assert(jj < iElevs[ii].size());
            float x = iElevs[ii][jj];
            float y = iGfield(ii,jj,e);
            float currLaf = iLafs[ii][jj];
            if(mLogTransform) {
               y = log(y);
            }
            if(Util::isValid(x) && Util::isValid(y) &&(!Util::isValid(currLaf) || currLaf >= mMinLafForElevGradient)) {
               meanXY += x*y;
               meanX  += x;
               meanY  += y;
               meanXX += x*x;
               counter++;
               // Found a new min
               if(!Util::isValid(min) || x < min)
                  min = x;
               // Found a new max
               if(!Util::isValid(max) || x > max)
                  max = x;
            }
         }
      }
      // Compute elevation difference within neighbourhood
      float elevDiff = Util::MV;
      if(Util::isValid(min) && Util::isValid(max)) {
         assert(max >= min);
         elevDiff = max - min;
      }

      // Use model gradient if:
      // 1) sufficient elevation difference in neighbourhood
      // 2) regression parameters are stable enough
      if(counter >= mMinNumPoints && Util::isValid(elevDiff) && elevDiff >= mMinElevDiff && meanXX != meanX*meanX) {
         // Estimate lapse rate
         meanXY /= counter;
         meanX  /= counter;
         meanY  /= counter;
         meanXX /= counter;
         if(meanXX - meanX*meanX != 0) {
            elevGradient = (meanXY - meanX*meanY)/(meanXX - meanX*meanX);
         }
      }
      else if(!mHasIssuedWarningUnstable) {
         std::stringstream ss;
         ss << "DownscalerGradient cannot compute gradient (unstable regression). Reverting to default gradient. Warning issued only once.";
         Util::warning(ss.str());
         mHasIssuedWarningUnstable = true;
      }
   }

   // Safety check
   if(!Util::isValid(elevGradient))
      elevGradient = mDefaultGradient;
   // Check against minimum and maximum gradients
   if(Util::isValid(mMinGradient) && elevGradient < mMinGradient)
      elevGradient = mMinGradient;
   if(Util::isValid(mMaxGradient) && elevGradient > mMaxGradient)
      elevGradient = mMaxGradient;
   return elevGradient;
}
#if 0
float DownscalerGradient2::get_laf_gradient(float iLaf, int iIcenter, int iJcenter, int e, vec2 iLafs, vec2 iElevs, const Field& iField, const File& iInput) const {
   // Compute the average temperature and LAF in a small radius
   float minLaf = Util::MV;
   float maxLaf = Util::MV;
   float minElev = Util::MV;
   float maxElev = Util::MV;
   float minT = Util::MV;
   float maxT = Util::MV;
   // std::vector<float> values(mSearchRadii.size(), Util::MV);
   float value = 0;
   int counter = 0;
   float totalWeight = 0;
   for(int r = 0; r < mLafSearchRadii.size(); r++) {
      int searchRadius = mLafSearchRadii[r];
      for(int ii = std::max(0, iIcenter-searchRadius); ii <= std::min(iInput.getNumLat()-1, iIcenter+searchRadius); ii++) {
         for(int jj = std::max(0, iJcenter-searchRadius); jj <= std::min(iInput.getNumLon()-1, iJcenter+searchRadius); jj++) {
            float currLaf  = iLafs[ii][jj];
            float currElev = iElevs[ii][jj];
            float currValue = iField(ii,jj,e);
            if(Util::isValid(currValue) && Util::isValid(currLaf)) {
               if(!Util::isValid(minLaf) || currLaf < minLaf) {
                  minLaf = currLaf;
                  minElev = currElev;
                  minT = currValue;
               }
               if(!Util::isValid(maxLaf) || currLaf > maxLaf) {
                  maxLaf = currLaf;
                  maxElev = currElev;
                  maxT = currValue;
               }
            }
         }
      }
      // std::cout << minLaf << " " << maxLaf << " " << minT << " " << maxT << std::endl;
      if(Util::isValid(maxLaf) && Util::isValid(minLaf) && (maxLaf - minLaf) >= mMinLafDiff) {
         float lafGradient = (maxT - minT) / (maxLaf - minLaf);
         // Altitude-adjust maxT so that it is on the same elevation as minT
         // Before computing the gradient
         if(Util::isValid(maxElev) && Util::isValid(minElev)) {
            float maxTcorr = maxT + mElevGradient * (maxElev - minElev);
            lafGradient = (maxTcorr - minT) / (maxLaf - minLaf);
         }
         float dLaf = currLaf - minLaf;
         // Safety check
         if(!Util::isValid(lafGradient))
            lafGradient = 0;
         // Check against minimum and maximum gradients
         if(Util::isValid(mMinGradient) && lafGradient < mMinGradient)
            lafGradient = mMinGradient;
         if(Util::isValid(mMaxGradient) && lafGradient > mMaxGradient)
            lafGradient = mMaxGradient;
         // std::cout << lafGradient << std::endl;
         // float currValue = nearestValue + dLaf * lafGradient;
         float currValue = minT + dLaf * lafGradient;

         // Altitude-adjust the current point
         if(Util::isValid(currElev) && Util::isValid(minElev)) {
            float dElev = currElev - minElev;
            currValue = currValue + dElev * mElevGradient;
         }
         // std::cout << currValue << std::endl;
         value += currValue * mWeights[r];
         totalWeight += mWeights[r];
         counter++;
         if(j == 43 && i == 56 && t == 0) {
            float dElev = currElev - minElev;
            float elevCorrection = dElev * mElevGradient;
            std::cout << minLaf << " " << maxLaf << " " << currLaf << "\n"
                      << minT << " " << maxT << " " << nearestValue << "\n"
                      << minElev << " " << maxElev << " " << currElev << "\n"
                      << currValue << " " << nearestValue << " " << value << "\n"
                      << lafGradient << " " << elevCorrection << std::endl;
         }
      }
   }
   // std::cout << value << std::endl;
   if(counter == 0 || !Util::isValid(value) || (Util::isValid(minAllowed) && value < minAllowed) || (Util::isValid(maxAllowed) && value > maxAllowed)) {
      value = nearestValue;
   }
   else {
      value = value / totalWeight;
   }

}

float DownscalerGradient2::calcElevationGradient() const {
   if(!Util::isValid(currElev) || !Util::isValid(nearestElev)) {
      // Can't adjust if we don't have an elevation, use nearest neighbour
      ofield(i,j,e) = ifield(Icenter,Jcenter,e);
   }
   else {
      float gradient = Util::MV;
      float totalValue = 0;
      float totalElev = 0;
      int counter = 0;
      int averagingRadius = 0;
      if(mAverageNeighbourhood)
         averagingRadius = mSearchRadius;
      for(int ii = std::max(0, Icenter-averagingRadius); ii <= std::min(iInput.getNumLat()-1, Icenter+averagingRadius); ii++) {
         for(int jj = std::max(0, Jcenter-averagingRadius); jj <= std::min(iInput.getNumLon()-1, Jcenter+averagingRadius); jj++) {
            float currValue = ifield(ii,jj,e);
            float currElev  = ielevs[ii][jj];
            if(Util::isValid(currValue) && Util::isValid(currElev)) {
               totalValue += currValue;
               totalElev  += currElev;
               counter++;
            }
         }
      }
      float baseValue = Util::MV;
      float baseElev  = Util::MV;
      if(counter > 0) {
         baseValue = totalValue / counter;
         baseElev  = totalElev  / counter;
      }
      float dElev = currElev - baseElev;

      if(Util::isValid(mConstantGradient)) {
         gradient = mConstantGradient;
      }
      else {
         /* Compute the model's gradient:
            The gradient is computed by using linear regression on forecast ~ elevation
            using all forecasts within a neighbourhood. To produce stable results, there
            is a requirement that the elevation within the neighbourhood has a large
            range (see mMinElevDiff).

            For bounded variables (e.g. wind speed), the gradient approach could cause
            forecasts to go outside its domain (e.g. negative winds). If this occurs,
            the nearest neighbour is used.
         */
         float meanXY  = 0; // elev*T
         float meanX   = 0; // elev
         float meanY   = 0; // T
         float meanXX  = 0; // elev*elev
         int   counter = 0;
         float min = Util::MV;
         float max = Util::MV;
         for(int ii = std::max(0, Icenter-mSearchRadius); ii <= std::min(iInput.getNumLat()-1, Icenter+mSearchRadius); ii++) {
            for(int jj = std::max(0, Jcenter-mSearchRadius); jj <= std::min(iInput.getNumLon()-1, Jcenter+mSearchRadius); jj++) {
               assert(ii < ielevs.size());
               assert(jj < ielevs[ii].size());
               float x = ielevs[ii][jj];
               float y = gfield(ii,jj,e);
               if(i == 772 && j == 659 && t == 0) {
                  std::cout << x << " " << y << std::endl;
               }
               if(mLogTransform) {
                  y = log(y);
               }
               if(Util::isValid(x) && Util::isValid(y)) {
                  meanXY += x*y;
                  meanX  += x;
                  meanY  += y;
                  meanXX += x*x;
                  counter++;
                  // Found a new min
                  if(!Util::isValid(min) || x < min)
                     min = x;
                  // Found a new max
                  if(!Util::isValid(max) || x > max)
                     max = x;
               }
            }
         }
         // Compute elevation difference within neighbourhood
         float elevDiff = Util::MV;
         if(Util::isValid(min) && Util::isValid(max)) {
            assert(max >= min);
            elevDiff = max - min;
         }

         // Use model gradient if:
         // 1) sufficient elevation difference in neighbourhood
         // 2) regression parameters are stable enough
         if(counter > 0 && Util::isValid(elevDiff) && elevDiff >= mMinElevDiff && meanXX != meanX*meanX) {
            // Estimate lapse rate
            meanXY /= counter;
            meanX  /= counter;
            meanY  /= counter;
            meanXX /= counter;
            if(meanXX - meanX*meanX != 0) {
               gradient = (meanXY - meanX*meanY)/(meanXX - meanX*meanX);
               if(i == 772 && j == 659 && t == 0) {
                  std::cout << olats[i][j] << " " << olons[i][j] << " " << counter << " " << meanX << " " << meanY << " " << meanXX << " " << gradient << std::endl;
               }
            }
         }
         else if(!mHasIssuedWarningUnstable) {
            std::stringstream ss;
            ss << "DownscalerGradient2 cannot compute gradient (unstable regression). Reverting to default gradient. Warning issued only once.";
            Util::warning(ss.str());
            mHasIssuedWarningUnstable = true;
         }
      }
      float value = Util::MV;
      if(mSaveGradient) {
         ofield(i, j, e) = gradient;
      }
      else {
         // Safety check
         if(!Util::isValid(gradient))
            gradien t= mDefaultGradient;
         // Check against minimum and maximum gradients
         if(Util::isValid(mMinGradient) && gradient < mMinGradient)
            gradient = mMinGradient;
         if(Util::isValid(mMaxGradient) && gradient > mMaxGradient)
            gradient = mMaxGradient;

         if(mLogTransform) {
            value = baseValue * exp(gradient * dElev);
         }
         else {
            value = baseValue + dElev * gradient;
         }
         if((Util::isValid(minAllowed) && value < minAllowed) || (Util::isValid(maxAllowed) && value > maxAllowed)) {
            // Use nearest neighbour if the gradient put us outside the bounds of the variable
            ofield(i,j,e) = baseValue;
         }
         else {
            ofield(i,j,e)  = value;
         }
      }
   }
}
#endif
float DownscalerGradient2::getConstantGradient() const {
   return mConstantGradient;
}
int DownscalerGradient2::getSearchRadius() const {
   return mSearchRadius;
}
float DownscalerGradient2::getMinElevDiff() const {
   return mMinElevDiff;
}
bool DownscalerGradient2::getLogTransform() const {
   return mLogTransform;
}
float DownscalerGradient2::getMinGradient() const {
   return mMinGradient;
}
float DownscalerGradient2::getMaxGradient() const {
   return mMaxGradient;
}
float DownscalerGradient2::getDefaultGradient() const {
   return mDefaultGradient;
}
std::string DownscalerGradient2::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-d gradient2", "Combines -d gradient and -d coastal.") << std::endl;
   ss << Util::formatDescription("   constantGradient=undef", "Fix gradient to this value. Units are units_of_variable / meter. Positive values mean an increase with height. If unspecified, computes the gradient by regression of points in a neighbourhood. The remaining options are ignored.") << std::endl;
   ss << Util::formatDescription("   searchRadius=3", "Compute gradient in a neighbourhood box of points within within +- radius in both east-west and north-south direction. Also use this neighbourhood to compute the vertical gradient.") << std::endl;
   ss << Util::formatDescription("   minElevDiff=30", "Minimum elevation range (in meters) required within the neighbourhood to compute gradient. Use nearest neighbour otherwise.") << std::endl;
   ss << Util::formatDescription("   lafRadius=3", "Compute gradient in a neighbourhood box of points within within +- radius in both east-west and north-south direction. Also use this neighbourhood to compute the vertical gradient.") << std::endl;
   ss << Util::formatDescription("   minLafDiff=30", "Minimum elevation range (in meters) required within the neighbourhood to compute gradient. Use nearest neighbour otherwise.") << std::endl;
   ss << Util::formatDescription("   logTransform=0", "Should the variable be log-transformed first? I.e should a linear gradient be applied to log(variable)? T = T(nn) * exp(gradient * dElev). Useful for pressure variables. Can be used together with constantGradient.") << std::endl;
   ss << Util::formatDescription("   defaultGradient=0", "If the gradient is not computable (too unstable), use this default gradient.") << std::endl;
   ss << Util::formatDescription("   minGradient=undef", "Do not allow gradient to be smaller than this value. If undefined, do not alter gradient.") << std::endl;
   ss << Util::formatDescription("   maxGradient=undef", "Do not allow gradient to be larger than this value. If undefined, do not alter gradient.") << std::endl;
   ss << Util::formatDescription("   averageNeighbourhood=0", "Should the average forecast and elevation within the search radius be used when determining what value to apply the gradient to?") << std::endl;
   ss << Util::formatDescription("   saveGradient=""", "Store the gradient instead of the value for debugging purposes. Use elev to store the elevation gradient; Use laf to store the laf gradient.") << std::endl;
   ss << Util::formatDescription("   minNumPoints=2", "Minimum number of points needed to compute an elevation gradient.") << std::endl;
   ss << Util::formatDescription("   minFracSeaPoints=0", "Minimum fraction of points in neighbourhood with a LAF < 0.5 to compute a LAF gradient.") << std::endl;
   ss << Util::formatDescription("   gradientVariable=T", "Which variable to use for the gradient?") << std::endl;
   return ss.str();
}

