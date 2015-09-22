#include "Gradient.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerGradient::DownscalerGradient(Variable::Type iVariable, const Options& iOptions) :
      Downscaler(iVariable, iOptions),
      mSearchRadius(3),
      mMinGradient(Util::MV),
      mMaxGradient(Util::MV),
      mLogTransform(false),
      mConstantGradient(Util::MV),
      mDefaultGradient(0),
      mMinElevDiff(30),
      mAverageNeighbourhood(false),
      mHasIssuedWarningUnstable(false) {

   iOptions.getValue("minGradient", mMinGradient);
   iOptions.getValue("maxGradient", mMaxGradient);
   iOptions.getValue("defaultGradient", mDefaultGradient);
   iOptions.getValue("searchRadius", mSearchRadius);
   iOptions.getValue("logTransform", mLogTransform);
   iOptions.getValue("constantGradient", mConstantGradient);
   iOptions.getValue("minElevDiff", mMinElevDiff);
   iOptions.getValue("averageNeighbourhood", mAverageNeighbourhood);
}

void DownscalerGradient::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumLat();
   int nLon = iOutput.getNumLon();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   vec2 ilats  = iInput.getLats();
   vec2 ilons  = iInput.getLons();
   vec2 ielevs = iInput.getElevs();
   vec2 olats  = iOutput.getLats();
   vec2 olons  = iOutput.getLons();
   vec2 oelevs = iOutput.getElevs();

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
               float currElev = oelevs[i][j];
               float nearestElev = ielevs[Icenter][Jcenter];
               if(!Util::isValid(currElev) || !Util::isValid(nearestElev)) {
                  // Can't adjust if we don't have an elevation, use nearest neighbour
                  ofield(i,j,e) = ifield(Icenter,Jcenter,e);
               }
               else {
                  float gradient = mDefaultGradient;
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
                           float y = ifield(ii,jj,e);
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
                        gradient = (meanXY - meanX*meanY)/(meanXX - meanX*meanX);
                     }
                     else if(!mHasIssuedWarningUnstable) {
                        std::stringstream ss;
                        ss << "DownscalerGradient cannot compute gradient (unstable regression). Reverting to default gradient. Warning issued only once.";
                        Util::warning(ss.str());
                        mHasIssuedWarningUnstable = true;
                     }
                     // Safety check
                     if(!Util::isValid(gradient))
                        gradient = 0;
                     // Check against minimum and maximum gradients
                     if(Util::isValid(mMinGradient) && gradient < mMinGradient)
                        gradient = mMinGradient;
                     if(Util::isValid(mMaxGradient) && gradient > mMaxGradient)
                        gradient = mMaxGradient;

                  }
                  float value = Util::MV;
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
      }
   }
}
float DownscalerGradient::getConstantGradient() const {
   return mConstantGradient;
}
int DownscalerGradient::getSearchRadius() const {
   return mSearchRadius;
}
float DownscalerGradient::getMinElevDiff() const {
   return mMinElevDiff;
}
bool DownscalerGradient::getLogTransform() const {
   return mLogTransform;
}
float DownscalerGradient::getMinGradient() const {
   return mMinGradient;
}
float DownscalerGradient::getMaxGradient() const {
   return mMaxGradient;
}
float DownscalerGradient::getDefaultGradient() const {
   return mDefaultGradient;
}
std::string DownscalerGradient::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-d gradient", "Adjusts the nearest neighbour based on the elevation difference to the output gridpoint and the gradient in the surrounding neighbourhood: T = T(nn) + gradient * dElev. If the gradient puts the forecast outside the domain of the variable (e.g. negative precipitation) then the nearest neighbour is used.") << std::endl;
   ss << Util::formatDescription("   constantGradient=undef", "Fix gradient to this value. Units are units_of_variable / meter. Positive values mean an increase with height. If unspecified, computes the gradient by regression of points in a neighbourhood. The remaining options are ignored.") << std::endl;
   ss << Util::formatDescription("   searchRadius=3", "Compute gradient in a neighbourhood box of points within within +- radius in both east-west and north-south direction. Also use this neighbourhood to compute the vertical gradient.") << std::endl;
   ss << Util::formatDescription("   minElevDiff=30", "Minimum elevation range (in meters) required within the neighbourhood to compute gradient. Use nearest neighbour otherwise.") << std::endl;
   ss << Util::formatDescription("   logTransform=0", "Should the variable be log-transformed first? I.e should a linear gradient be applied to log(variable)? T = T(nn) * exp(gradient * dElev). Useful for pressure variables. Can be used together with constantGradient.") << std::endl;
   ss << Util::formatDescription("   defaultGradient=0", "If the gradient is not computable (too unstable), use this default gradient.") << std::endl;
   ss << Util::formatDescription("   minGradient=undef", "Do not allow gradient to be smaller than this value. If undefined, do not alter gradient.") << std::endl;
   ss << Util::formatDescription("   maxGradient=undef", "Do not allow gradient to be larger than this value. If undefined, do not alter gradient.") << std::endl;
   ss << Util::formatDescription("   averageNeighbourhood=0", "Should the average forecast and elevation within the search radius be used when determining what value to apply the gradient to?") << std::endl;
   return ss.str();
}

