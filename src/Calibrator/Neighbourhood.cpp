#include "Neighbourhood.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
CalibratorNeighbourhood::CalibratorNeighbourhood(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mRadius(3),
      mStatType(Util::StatTypeMean),
      mFast(true),
      mQuantile(Util::MV) {
   iOptions.getValue("radius", mRadius);
   iOptions.getValue("fast", mFast);
   if(mRadius < 0) {
      std::stringstream ss;
      ss << "CalibratorNeighbourhood: 'radius' (" << mRadius << ") must be >= 0";
      Util::error(ss.str());
   }

   std::string op;
   if(iOptions.getValue("stat", op)) {
      bool status = Util::getStatType(op, mStatType);
      if(!status) {
         std::stringstream ss;
         ss << "Could not recognize stat=" << op;
         Util::error(ss.str());
      }
   }
   if(mStatType == Util::StatTypeQuantile) {
      iOptions.getRequiredValue("quantile", mQuantile);
      if(!Util::isValid(mQuantile) || mQuantile < 0 || mQuantile > 1) {
         Util::error("'quantile' must be on the interval [0,1]");
      }
   }
   iOptions.check();
}

bool CalibratorNeighbourhood::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& precip = *iFile.getField(mVariable, t);
      Field precipRaw = precip;

      int radius = mRadius;
      if(iParameterFile != NULL) {
         if(iParameterFile->isLocationDependent()) {
            Util::error("Cannot use a location dependent parameter file for CalibratorNeighbourhood");
         }
         radius = iParameterFile->getParameters(t)[0];
      }

      if(mFast && (mStatType == Util::StatTypeMean || mStatType == Util::StatTypeSum)) {
         for(int e = 0; e < nEns; e++) {
            vec2 values;
            vec2 counts;
            values.resize(nLat);
            counts.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               values[i].resize(nLon, 0);
               counts[i].resize(nLon, 0);
            }
            // Compute accumulated values
            for(int i = 0; i < nLat; i++) {
               for(int j = 0; j < nLon; j++) {
                  float value = precipRaw(i, j, e);
                  if(j == 0 && i == 0) {
                     // Lower corner
                     if(Util::isValid(value)) {
                        values[i][j] = precipRaw(i, j, e);
                        counts[i][j] = 1;
                     }
                  }
                  else if(j == 0) {
                     // Lower row
                     if(Util::isValid(value)) {
                        values[i][j] = values[i-1][j] + precipRaw(i,j,e);
                        counts[i][j] = counts[i-1][j] + 1;
                     }
                     else {
                        values[i][j] = values[i-1][j];
                        counts[i][j] = counts[i-1][j];
                     }
                  }
                  else if(i == 0) {
                     // Left column
                     if(Util::isValid(value)) {
                        values[i][j] = values[i][j-1] + precipRaw(i,j,e);
                        counts[i][j] = counts[i][j-1] + 1;
                     }
                     else {
                        values[i][j] = values[i][j-1];
                        counts[i][j] = counts[i][j-1];
                     }

                  }
                  else {
                     if(Util::isValid(value)) {
                        values[i][j] = values[i][j-1] + values[i-1][j] - values[i-1][j-1] + precipRaw(i,j,e);
                        counts[i][j] = counts[i][j-1] + counts[i-1][j] - counts[i-1][j-1] + 1;
                     }
                     else {
                        values[i][j] = values[i][j-1] + values[i-1][j] - values[i-1][j-1];
                        counts[i][j] = counts[i][j-1] + counts[i-1][j] - counts[i-1][j-1];
                     }
                  }
               }
            }
            // Put neighbourhood into vector
            #pragma omp parallel for
            for(int i = 0; i < nLat; i++) {
               for(int j = 0; j < nLon; j++) {
                  int i1 = std::min(nLat-1, i + mRadius);
                  int j1 = std::min(nLon-1, j + mRadius);
                  int i0 = i - mRadius - 1;
                  int j0 = j - mRadius - 1;
                  float value11 = values[i1][j1];
                  float value00 = 0;
                  float value10 = 0;
                  float value01 = 0;
                  int count11 = counts[i1][j1];
                  int count00 = 0;
                  int count10 = 0;
                  int count01 = 0;
                  if(i0 >= 0 && j0 >= 0) {
                     value00 = values[i0][j0];
                     value10 = values[i1][j0];
                     value01 = values[i0][j1];
                     count00 = counts[i0][j0];
                     count10 = counts[i1][j0];
                     count01 = counts[i0][j1];
                  }
                  else if(j0 >= 0) {
                     value10 = values[i1][j0];
                     count10 = counts[i1][j0];
                  }
                  else if(i0 >= 0) {
                     value01 = values[i0][j1];
                     count01 = counts[i0][j1];
                  }
                  float value = value11 + value00 - value10 - value01;
                  int count = count11 + count00 - count10 - count01;
                  if(count > 0) {
                     if(mStatType == Util::StatTypeMean) {
                        value /= count;
                     }
                     precip(i,j,e) = value;
                  }
               }
            }
         }
      }
      else {
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  // Put neighbourhood into vector
                  std::vector<float> neighbourhood;
                  int Ni = std::min(nLat-1, i+radius) - std::max(0, i-radius) + 1;
                  int Nj = std::min(nLon-1, j+radius) - std::max(0, j-radius) + 1;
                  assert(Ni > 0);
                  assert(Nj > 0);
                  neighbourhood.resize(Ni*Nj, Util::MV);
                  int index = 0;
                  for(int ii = std::max(0, i-radius); ii <= std::min(nLat-1, i+radius); ii++) {
                     for(int jj = std::max(0, j-radius); jj <= std::min(nLon-1, j+radius); jj++) {
                        float value = precipRaw(ii,jj,e);
                        assert(index < Ni*Nj);
                        neighbourhood[index] = value;
                        index++;
                     }
                  }
                  assert(index == Ni*Nj);
                  precip(i,j,e) = Util::calculateStat(neighbourhood, mStatType, mQuantile);
               }
            }
         }
      }
   }
   return true;
}

std::string CalibratorNeighbourhood::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c neighbourhood", "Applies a statistical operator on a neighbourhood (example by averaging across a neighbourhood thereby smoothing the field).") << std::endl;
      ss << Util::formatDescription("   radius=3", "Use gridpoints within this number of points within in both east-west and north-south direction. The radius can alternatively be specified using a location-independent parameter file, with one parameter.") << std::endl;
      ss << Util::formatDescription("   stat=mean", "What statistical operator should be applied to the neighbourhood? One of 'mean', 'median', 'min', 'max', 'quantile', 'std', or 'sum'. 'std' is the population standard deviation.") << std::endl;
      ss << Util::formatDescription("   quantile=undef", "If stat=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
      ss << Util::formatDescription("   fast=1", "Use shortcuts to compute 'mean' or 'sum faster.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c neighbourhood", "Applies a statistical operator on a neighbourhood") << std::endl;
   return ss.str();
}
