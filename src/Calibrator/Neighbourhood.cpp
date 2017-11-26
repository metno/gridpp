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
      mQuantile(Util::MV) {
   iOptions.getValue("radius", mRadius);
   if(mRadius < 0) {
      std::stringstream ss;
      ss << "CalibratorNeighbourhood: 'radius' (" << mRadius << ") must be >= 0";
      Util::error(ss.str());
   }

   std::string op;
   if(iOptions.getValue("stat", op)) {
      if(op == "mean") {
         mStatType = Util::StatTypeMean;
      }
      else if(op == "min") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 0;
      }
      else if(op == "max") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 1;
      }
      else if(op == "median") {
         mStatType = Util::StatTypeQuantile;
         mQuantile = 0.5;
      }
      else if(op == "std") {
         mStatType = Util::StatTypeStd;
      }
      else if(op == "quantile"){
         mStatType = Util::StatTypeQuantile;
         if(!iOptions.getValue("quantile", mQuantile)) {
            Util::error("CalibratorNeighbourhood: option 'quantile' is required");
         }
         if(!Util::isValid(mQuantile) || mQuantile < 0 || mQuantile > 1) {
            Util::error("CalibratorNeighbourhood: 'quantile' must be on the interval [0,1]");
         }
      }
      else {
         Util::error("CalibratorNeighbourhood: Unrecognized value for 'stat'");
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
   return true;
}

std::string CalibratorNeighbourhood::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c neighbourhood", "Applies a statistical operator on a neighbourhood (example by averaging across a neighbourhood thereby smoothing the field).") << std::endl;
   ss << Util::formatDescription("   radius=3", "Use gridpoints within this number of points within in both east-west and north-south direction. The radius can alternatively be specified using a location-independent parameter file, with one parameter.") << std::endl;
   ss << Util::formatDescription("   stat=mean", "What statistical operator should be applied to the neighbourhood? One of 'mean', 'median', 'min', 'max', 'std', or 'quantile'. 'std' is the population standard deviation.") << std::endl;
   ss << Util::formatDescription("   quantile=undef", "If stat=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
   return ss.str();
}
