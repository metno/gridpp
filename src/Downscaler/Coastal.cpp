#include "Coastal.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerCoastal::DownscalerCoastal(Variable::Type iVariable, const Options& iOptions) :
      Downscaler(iVariable, iOptions),
      mSearchRadius(3),
      mMinGradient(Util::MV),
      mMaxGradient(Util::MV),
      mMinLafDiff(0.1) {

   iOptions.getValue("minGradient", mMinGradient);
   iOptions.getValue("maxGradient", mMaxGradient);
   iOptions.getValue("searchRadius", mSearchRadius);
   iOptions.getValue("minLafDiff", mMinLafDiff);
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
               float nearestLaf = ilafs[Icenter][Jcenter];
               if(!Util::isValid(currLaf) || !Util::isValid(nearestLaf)) {
                  // Can't adjust if we don't have LAF, use nearest neighbour
                  ofield(i,j,e) = ifield(Icenter,Jcenter,e);
               }
               else {
                  float gradient = Util::MV;
                  float totalValue = 0;
                  float totalLaf = 0;
                  int counter = 0;
                  int averagingRadius = 0;
                  /* Here is the formula:
                     T = Tmin + (Tmax  - Tmin) / (LAFmax - LAFmin) * (LAF - Tmin)
                  */
                  // Compute the average temperature and LAF in a small radius
                  float lafMin = Util::MV;
                  float lafMax = Util::MV;
                  float tempMin = Util::MV;
                  float tempMax = Util::MV;
                  for(int ii = std::max(0, Icenter-averagingRadius); ii <= std::min(iInput.getNumLat()-1, Icenter+averagingRadius); ii++) {
                     for(int jj = std::max(0, Jcenter-averagingRadius); jj <= std::min(iInput.getNumLon()-1, Jcenter+averagingRadius); jj++) {
                        float currValue = ifield(ii,jj,e);
                        float currLaf  = ilafs[ii][jj];
                        if(Util::isValid(currValue) && Util::isValid(currLaf)) {
                           totalValue += currValue;
                           totalLaf  += currLaf;
                           counter++;
                        }
                     }
                  }
                  float baseValue = Util::MV;
                  float baseLaf  = Util::MV;
                  if(counter > 0) {
                     baseValue = totalValue / counter;
                     baseLaf  = totalLaf  / counter;
                  }
                  float dLaf = currLaf - baseLaf;

                  /* Compute the model's gradient:
                     The gradient is computed by using linear regression on forecast ~ laf
                     using all forecasts within a neighbourhood. To produce stable results, there
                     is a requirement that the elevation within the neighbourhood has a large
                     range (see mMinLafDiff).

                     For bounded variables (e.g. wind speed), the gradient approach could cause
                     forecasts to go outside its domain (e.g. negative winds). If this occurs,
                     the nearest neighbour is used.
                  */
                  float meanXY  = 0; // elev*T
                  float meanX   = 0; // elev
                  float meanY   = 0; // T
                  float meanXX  = 0; // elev*elev

                  counter = 0;
                  float min = Util::MV;
                  float max = Util::MV;
                  for(int ii = std::max(0, Icenter-mSearchRadius); ii <= std::min(iInput.getNumLat()-1, Icenter+mSearchRadius); ii++) {
                     for(int jj = std::max(0, Jcenter-mSearchRadius); jj <= std::min(iInput.getNumLon()-1, Jcenter+mSearchRadius); jj++) {
                        assert(ii < ilafs.size());
                        assert(jj < ilafs[ii].size());
                        float x = ilafs[ii][jj];
                        float y = ifield(ii,jj,e);
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
                  // Compute laf difference within neighbourhood
                  float lafDiff = Util::MV;
                  if(Util::isValid(min) && Util::isValid(max)) {
                     assert(max >= min);
                     lafDiff = max - min;
                  }

                  // Use model gradient if:
                  // 1) sufficient laf difference in neighbourhood
                  // 2) regression parameters are stable enough
                  if(counter > 0 && Util::isValid(lafDiff) && lafDiff >= mMinLafDiff && meanXX != meanX*meanX) {
                     // Estimate lapse rate
                     meanXY /= counter;
                     meanX  /= counter;
                     meanY  /= counter;
                     meanXX /= counter;
                     gradient = (meanXY - meanX*meanY)/(meanXX - meanX*meanX);
                  }
                  // Safety check
                  if(!Util::isValid(gradient))
                     gradient = 0;
                  // Check against minimum and maximum gradients
                  if(Util::isValid(mMinGradient) && gradient < mMinGradient)
                     gradient = mMinGradient;
                  if(Util::isValid(mMaxGradient) && gradient > mMaxGradient)
                     gradient = mMaxGradient;
                  float value = Util::MV;
                  value = baseValue + dLaf * gradient;
                  if((Util::isValid(minAllowed) && value < minAllowed) || (Util::isValid(maxAllowed) && value > maxAllowed)) {
                     // Use nearest neighbour if the gradient put us outside the bounds of the variable
                     ofield(i,j,e) = baseValue;
                  }
                  else {
                     // ofield(i,j,e)  = gradient + 273.15;
                     ofield(i,j,e)  = value;
                  }
               }
            }
         }
      }
   }
}
int DownscalerCoastal::getSearchRadius() const {
   return mSearchRadius;
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

