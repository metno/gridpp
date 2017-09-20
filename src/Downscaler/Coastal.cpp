#include "Coastal.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerCoastal::DownscalerCoastal(Variable::Type iVariable, const Options& iOptions) :
      Downscaler(iVariable, iOptions),
      mMinGradient(Util::MV),
      mMaxGradient(Util::MV),
      mLafRadius(1),
      mElevGradient(-0.0065),
      mMinLafDiff(0.1) {

   if(!iOptions.getValues("searchRadii", mSearchRadii)) {
      mSearchRadii.push_back(1);
      mSearchRadii.push_back(2);
      mSearchRadii.push_back(3);
   }

   if(!iOptions.getValues("weights", mWeights)) {
      mWeights.push_back(0.5);
      mWeights.push_back(0.3);
      mWeights.push_back(0.2);
   }

   assert(mSearchRadii.size() == mWeights.size());

   iOptions.getValue("minGradient", mMinGradient);
   iOptions.getValue("maxGradient", mMaxGradient);
   iOptions.getValue("minLafDiff", mMinLafDiff);
   iOptions.getValue("lafRadius", mLafRadius);
   iOptions.getValue("elevGradient", mElevGradient);
}

void DownscalerCoastal::downscaleCore(const File& iInput, File& iOutput) const {
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

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mVariable, t);
      Field& ofield = *iOutput.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int Icenter = nearestI[i][j];
            int Jcenter = nearestJ[i][j];
            assert(Icenter < ielevs.size());
            assert(Jcenter < ielevs[Icenter].size());
            for(int e = 0; e < nEns; e++) {
               float currLaf = olafs[i][j];
               float currElev = oelevs[i][j];
               if(mLafRadius > 0) {
                  float total = 0;
                  int counter = 0;
                  for(int ii = std::max(0, i-mLafRadius); ii <= std::min(nLat-1, i+mLafRadius); ii++) {
                     for(int jj = std::max(0, j-mLafRadius); jj <= std::min(nLon-1, j+mLafRadius); jj++) {
                        if(Util::isValid(olafs[ii][jj])) {
                           total += olafs[ii][jj];
                           counter++;
                        }
                     }
                  }
                  if(counter > 0)
                     currLaf = total / counter;
               }
               float nearestLaf = ilafs[Icenter][Jcenter];
               float nearestValue = ifield(Icenter,Jcenter,e);
               if(!Util::isValid(currLaf) || !Util::isValid(nearestLaf)) {
                  // Can't adjust if we don't have LAF, use nearest neighbour
                  ofield(i,j,e) = nearestValue;
               }
               else {
                  /* Here is the formula:
                     value = minT + (maxT  - minT) / (maxLaf - minLaf) * (laf - minLaf)
                     value = minT + gradient * dLaf
                  */
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
                  for(int r = 0; r < mSearchRadii.size(); r++) {
                     int searchRadius = mSearchRadii[r];
                     for(int ii = std::max(0, Icenter-searchRadius); ii <= std::min(iInput.getNumLat()-1, Icenter+searchRadius); ii++) {
                        for(int jj = std::max(0, Jcenter-searchRadius); jj <= std::min(iInput.getNumLon()-1, Jcenter+searchRadius); jj++) {
                           float currLaf  = ilafs[ii][jj];
                           float currElev = ielevs[ii][jj];
                           float currValue = ifield(ii,jj,e);
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
                  ofield(i,j,e)  = value;
               }
            }
         }
      }
   }
}
std::vector<int> DownscalerCoastal::getSearchRadii() const {
   return mSearchRadii;
}
float DownscalerCoastal::getMinLafDiff() const {
   return mMinLafDiff;
}
float DownscalerCoastal::getMinGradient() const {
   return mMinGradient;
}
float DownscalerCoastal::getMaxGradient() const {
   return mMaxGradient;
}
std::string DownscalerCoastal::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-d coastal", "Computes the gradient with respect to land area fraction and interpolates based on this.") << std::endl;
   ss << Util::formatDescription("   constantGradient=undef", "Fix gradient to this value. Units are units_of_variable / meter. Positive values mean an increase with height. If unspecified, computes the gradient by regression of points in a neighbourhood. The remaining options are ignored.") << std::endl;
   ss << Util::formatDescription("   searchRadius=3", "Compute gradient in a neighbourhood box of points within within +- radius in both east-west and north-south direction. Also use this neighbourhood to compute the vertical gradient.") << std::endl;
   ss << Util::formatDescription("   minLafDiff=30", "Minimum elevation range (in meters) required within the neighbourhood to compute gradient. Use nearest neighbour otherwise.") << std::endl;
   ss << Util::formatDescription("   logTransform=0", "Should the variable be log-transformed first? I.e should a linear gradient be applied to log(variable)? T = T(nn) * exp(gradient * dElev). Useful for pressure variables. Can be used together with constantGradient.") << std::endl;
   ss << Util::formatDescription("   defaultGradient=0", "If the gradient is not computable (too unstable), use this default gradient.") << std::endl;
   ss << Util::formatDescription("   minGradient=undef", "Do not allow gradient to be smaller than this value. If undefined, do not alter gradient.") << std::endl;
   ss << Util::formatDescription("   maxGradient=undef", "Do not allow gradient to be larger than this value. If undefined, do not alter gradient.") << std::endl;
   ss << Util::formatDescription("   averageNeighbourhood=0", "Should the average forecast and elevation within the search radius be used when determining what value to apply the gradient to?") << std::endl;
   return ss.str();
}

