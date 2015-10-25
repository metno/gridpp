#include "KalmanFilter.h"
#include "Util.h"
#include "Parameters.h"
#include "Options.h"
#include "Location.h"
#include "Field.h"
#include "ParameterFile/ParameterFile.h"
#include "Downscaler/Downscaler.h"
#include "File/File.h"
#include <math.h>
#include <assert.h>

KalmanFilter::KalmanFilter(Variable::Type iVariable, const Options& iOptions) :
      mVariable(iVariable),
      mHourlyCorr(0.93),
      mV(2.0),
      mRatio(0.06),
      mDim(8) {
   iOptions.getValue("hourlyCorr",  mHourlyCorr);
   iOptions.getValue("v",  mV);
   iOptions.getValue("ratio",  mRatio);
   iOptions.getValue("dim",  mDim);

   if(mHourlyCorr < 0 || mHourlyCorr > 1)
      Util::error("'hourlyCorr' must be between 0 and 1 inclusive");
   if(mV < 0)
      Util::error("'v' must be >= 0");
   if(mRatio < 0 || mRatio > 1)
      Util::error("'ratio' must be between 0 and 1 inclusive");
   if(mDim < 1 || mDim > 24)
      Util::error("'dim' must be between 1 and 24 inclusive");
}

bool KalmanFilter::writeBiasFile(const File& iFcst,
      const File& iObs,
      int iCurrDate,
      int iStartTime,
      int iEndTime,
      const ParameterFile* iDbIn,
      ParameterFile* iDbOut,
      ParameterFile* iBiasFile) {
   vec2 flats = iFcst.getLats();
   vec2 flons = iFcst.getLons();
   vec2 felevs = iFcst.getElevs();
   vec2 olats = iObs.getLats();
   vec2 olons = iObs.getLons();
   vec2 oelevs = iObs.getElevs();

   // Which forecast index are we using to update?
   std::vector<std::vector<float> > W = getW();

   vec2Int I,J;
   Downscaler::getNearestNeighbour(iFcst, iObs, I, J);

   // Loop over locations
   for(int oi = 0; oi < iObs.getNumLat(); oi++) {
      for(int oj = 0; oj < iObs.getNumLon(); oj++) {
         Location obsLoc(olats[oi][oj], olons[oi][oj], oelevs[oi][oj]);
         int Inearest = I[oi][oj];
         int Jnearest = J[oi][oj];
         bool inDb = true;
         if(inDb) {
            int iI, iJ;
            KalmanParameters par = initialize(); // Initialize to empty
            if(iDbIn != NULL) {
               const Parameters& rawPar = iDbIn->getParameters(0, obsLoc); // Same parameters for all hours
               par = KalmanParameters(rawPar);
            }
            std::vector<float> x = par.x;
            std::vector<float> k = par.k;
            std::vector<std::vector<float> > P= par.P;

            // Iterate forward
            for(int t = iStartTime; t <= iEndTime; t++) {
               std::cout << "Iteration: " << t << std::endl;
               int ftime = 0; // Forecast time
               int findex = t;
               float fcst = (*iFcst.getField(mVariable, findex))(Inearest,Jnearest,0);
               float obs = (*iObs.getField(mVariable, findex))(oi, oj,0);

               // Compute Pt|t-1. This increases the covariance.
               for(int i = 0; i < mDim; i++) {
                  for(int j = 0; j < mDim; j++) {
                     P[i][j] += W[i][j];
                  }
               }
               if(Util::isValid(fcst) && Util::isValid(obs)) {
                  float y = obs - fcst;
                  int dindex = getParameterIndex(t);
                  assert(dindex < k.size());
                  assert(dindex < P.size());

                  // Compute Kt
                  for(int i = 0; i < mDim; i++) {
                     k[i] = P[dindex][i] / (P[dindex][dindex] + mV);
                  }
                  // Compute Pt|t
                  for(int i = 0; i < mDim; i++) {
                     for(int j = 0; j < mDim; j++) {
                        P[i][j] = (1 - k[dindex])*P[i][j]; // ? i or j for k?
                     }
                  }
                  // Compute xt|t
                  float x_dindex = x[dindex];
                  for(int i = 0; i < mDim; i++) {
                     x[i] = x[i] + k[i]*(y - x_dindex);
                  }
                  std::cout << "Iteration " << t << " " << dindex << " " << y;
                  for(int i = 0; i < mDim; i++) {
                     std::cout << " " << x[i];
                  }
                  std::cout << std::endl;
                  std::cout << "# " << t << " " << x[0] << " " << k[0] << " " << P[0][0] << " " << P[1][0] << " " << y << std::endl;
               }
               else {
                  // Missing obs or forecast. P has already been increased, do not update x.
                  // TODO: Does the kalman gain need to be updated?
               }
            }
            
            // Store the new parameters
            KalmanParameters parNew = par;
            parNew.x = x;
            parNew.k = k;
            parNew.P = P;

            if(iDbOut != NULL) {
               iDbOut->setParameters(parNew.toParameters(), 0, obsLoc);
            }
            if(iBiasFile != NULL) {
               std::vector<float> values(1,0);
               int H = iFcst.getNumTime(); // TODO: Figure out which hours to write which parameters to
               for(int h = 0; h < H; h++) {
                  int index = round((float) h*mDim/24);
                  values[0] = parNew.x[index];
                  iBiasFile->setParameters(Parameters(values), h, obsLoc);
               }
            }
            std::cout << "gain bias" << std::endl;
            for(int i = 0; i < mDim; i++) {
               std::cout << k[i] << " " << parNew.x[i];
               for(int j = i; j < i+1; j++) {
                  std::cout << " " << P[i][j];
               }
               std::cout << std::endl;
            }

         }
         else {

         }
      }
   }
   if(iDbOut != NULL) {
      iDbOut->write();
   }
   if(iBiasFile != NULL) {
      iBiasFile->write();
   }
}

std::vector<std::vector<float> > KalmanFilter::getW() const {
   std::vector<std::vector<float> > vec;
   vec.resize(mDim);
   int hoursBetween = 24 / mDim;
   float factor = 100;
   if(mV != 0)
      factor = mRatio * mV;
   for(int i = 0; i < mDim; i++) {
      vec[i].resize(mDim, 0);
      for(int j = 0; j < mDim; j++) {
         int diff = std::min(abs(i-j), mDim - abs(i-j));
         float weight = pow(mHourlyCorr, (float) diff*hoursBetween);
         vec[i][j] = weight * factor;
      }
      // assert(vec[i][i] ==  1);
   }
   return vec;
}

KalmanParameters::KalmanParameters(const Parameters& iParameters) {
   int size = iParameters.size();
   int N = sqrt(1+size) - 1;
   int index = 0;
   for(int i = 0; i < N; i++) {
      x.push_back(iParameters[i]);
      index++;
   }
   for(int i = 0; i < N; i++) {
      k.push_back(iParameters[N+i]);
      index++;
   }
   P.resize(N);
   for(int i = 0; i < N; i++) {
      P[i].resize(N);
      for(int j = 0; j < N; j++) {
         P[i][j] = iParameters[index];
         index++;
      }
   }
}

/*
KalmanParameters::KalmanParameters(int N) {
   int index = 0;
   x.resize(N, 0);
   k.resize(N, 1);
   P.resize(N);
   for(int i = 0; i < N; i++) {
      P[i].resize(N,1000);
   }
}
*/
KalmanParameters KalmanFilter::initialize() const {
   int index = 0;
   std::vector<float> x;
   std::vector<float> k;
   std::vector<std::vector<float> > P = getW();
   for(int i = 0; i < P.size(); i++) {
      for(int j = 0; j < P.size(); j++) {
         // How should this be initialized?
         P[i][j] = P[i][j] / mRatio;
      }
   }
   x.resize(mDim, 0);
   k.resize(mDim, 0);
   return KalmanParameters(x, k, P);
}

KalmanParameters::KalmanParameters(std::vector<float> x, std::vector<float> k, std::vector<std::vector<float> > P) {
   this->x = x;
   this->k = k;
   this->P = P;
}

Parameters KalmanParameters::toParameters() const {
   std::vector<float> values;
   for(int i = 0; i < x.size(); i++) {
      values.push_back(x[i]);
   }
   for(int i = 0; i < k.size(); i++) {
      values.push_back(k[i]);
   }
   for(int i = 0; i < P.size(); i++) {
      for(int j = 0; j < P.size(); j++) {
         values.push_back(P[i][j]);
      }
   }
   Parameters par(values);
   return par;
}

int KalmanFilter::getParameterIndex(int iTime) const {
   int index = round((float) iTime * mDim / 24);
   index = index % mDim;
   assert(index >= 0 && index < mDim);
   return index;
}

std::string KalmanFilter::description() {
   std::stringstream ss;
   ss << Util::formatDescription("   hourlyCorr=0.93", "What is the correlation of bias between two neighbouring hours?") << std::endl;
   ss << Util::formatDescription("   v=2", "What is the standard deviation of the bias?") << std::endl;
   ss << Util::formatDescription("   ratio=0.06", "What is the ratio of the standard deviation of the change in the true bias to the standard deviation of the change of the measured bias? Higher values lead to a faster changing filter.") << std::endl;
   return ss.str();
}
