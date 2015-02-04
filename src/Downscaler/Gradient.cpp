#include "Gradient.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerGradient::DownscalerGradient(Variable::Type iVariable) :
      Downscaler(iVariable),
      mSearchRadius(3),
      mMinGradient(-10),
      mMaxGradient(10),
      mConstGradient(Util::MV),
      mMinElevDiff(30) {
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
   getNearestNeighbourFast(iInput, iOutput, nearestI, nearestJ);

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
                  float dElev = currElev - nearestElev;
                  float gradient = 0;
                  if(!Util::isValid(mConstGradient)) {
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

                     if(counter > 0 && Util::isValid(elevDiff) && elevDiff >= mMinElevDiff && meanXX != meanX*meanX) {
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
                  }
                  else {
                     gradient = mConstGradient;
                  }
                  float value = ifield(Icenter,Jcenter,e) + dElev * gradient;
                  if((Util::isValid(minAllowed) && value < minAllowed) || (Util::isValid(maxAllowed) && value > maxAllowed)) {
                     // Use nearest neighbour if we the gradient put us outside the bounds of the variable
                     ofield(i,j,e) = (ifield)(Icenter, Jcenter, e);
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
void DownscalerGradient::setConstantGradient(float iGradient) {
   if(!Util::isValid(iGradient)) {
      std::stringstream ss;
      ss << "DownscalerGradient: constant gradient must be a valid number";
      Util::error(ss.str());
   }
   mConstGradient = iGradient;
}
float DownscalerGradient::getConstantGradient() const {
   return mConstGradient;
}
void DownscalerGradient::setSearchRadius(int iNumPoints) {
   if(!Util::isValid(iNumPoints) || iNumPoints <= 0) {
      std::stringstream ss;
      ss << "DownscalerGradient: search radius must be >= 1";
      Util::error(ss.str());
   }
   mSearchRadius = iNumPoints;
}

int DownscalerGradient::getSearchRadius() const {
   return mSearchRadius;
}
void DownscalerGradient::setMinElevDiff(float iMinElevDiff) {
   mMinElevDiff = iMinElevDiff;
}
float DownscalerGradient::getMinElevDiff() const {
   return mMinElevDiff;
}

std::string DownscalerGradient::description() {
   std::stringstream ss;
   ss << "   -d gradient                  Adjusts the nearest neighbour based on the elevation difference" << std::endl;
   ss << "                                to the output gridpoint. If the gradient puts the forecast outside" << std::endl;
   ss << "                                to domain of the variable (e.g. negative precipitation) then the " << std::endl;
   ss << "                                nearest neighbour is used." << std::endl;
   ss << "      constantGradient=undef    Fix gradient to this value. If unspecified, computes the gradient" << std::endl;
   ss << "                                by linear regression of points in a neighbourhood." << std::endl;
   ss << "      searchRadius=3            Compute gradient in a neighbourhood box of points within +- radius" << std::endl;
   ss << "                                in both east-west and north-south direction." << std::endl;
   ss << "      minElevDiff=30            Minimum elevation range (in meters) required within the neighbourhood" << std::endl;
   ss << "                                to compute gradient. Use nearest neighbour otherwise." << std::endl;
   return ss.str();
}
