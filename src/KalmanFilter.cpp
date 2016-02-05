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
      mPinit(0.5),
      mElevGradient(-0.0065),
      mDim(8) {
   iOptions.getValue("hourlyCorr",  mHourlyCorr);
   iOptions.getValue("v",  mV);
   iOptions.getValue("ratio",  mRatio);
   iOptions.getValue("elevGradient",  mElevGradient);
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
      int iTimeStep,
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
   vec2Int I,J;
   Downscaler::getNearestNeighbour(iFcst, iObs, I, J);

   // Loop over locations
   for(int oi = 0; oi < iObs.getNumLat(); oi++) {
      for(int oj = 0; oj < iObs.getNumLon(); oj++) {
         Location obsLoc(olats[oi][oj], olons[oi][oj], oelevs[oi][oj]);
         int Inearest = I[oi][oj];
         int Jnearest = J[oi][oj];
         float oelev = oelevs[oi][oj];
         float felev = felevs[Inearest][Jnearest];

         // Compute the most recent bias
         float elevCorrection = mElevGradient * (oelev - felev);
         float fcst = (*iFcst.getField(mVariable, iTimeStep))(Inearest,Jnearest,0);
         fcst += elevCorrection;
         float obs = (*iObs.getField(mVariable, iTimeStep))(oi, oj,0);
         float bias = Util::MV;
         if(Util::isValid(obs) && Util::isValid(fcst))
            bias = obs - fcst;

         // Get KF parameters
         KalmanParameters par = initialize(); // Initialize to empty
         if(iDbIn != NULL) {
            const Parameters& rawPar = iDbIn->getParameters(0, obsLoc); // Same parameters for all hours
            par = KalmanParameters(rawPar);
         }

         // Update forecasts
         KalmanParameters parNew = update(bias, iTimeStep, par);

         // Store parameters and bias
         if(iDbOut != NULL) {
            iDbOut->setParameters(parNew.toParameters(), 0, obsLoc);
         }
         if(iBiasFile != NULL) {
            std::vector<float> values(1,0);
            // Figure out which hours to write which parameters to
            std::vector<double> times = iFcst.getTimes();
            double refTime = iFcst.getReferenceTime();
            for(int h = 0; h < times.size(); h++) {
               int sec = times[h] - refTime;
               int hour = sec / 3600 % 24;
               hour = h;
               int index = (int) round((float) hour*mDim/24) % mDim;
               values[0] = parNew.x[index];
               assert(parNew.x.size() > index);
               iBiasFile->setParameters(Parameters(values), h, obsLoc);
            }
         }
      }
   }
   if(iDbOut != NULL) {
      iDbOut->write();
   }
   if(iBiasFile != NULL) {
      iBiasFile->write();
   }
   return true;
}

vec2 KalmanFilter::getW() const {
   // Compute if it does not exist
   if(mW.size() == 0) {
      mW.resize(mDim);
      int hoursBetween = 24 / mDim;
      float factor = 100;
      if(mV != 0)
         factor = mRatio * mV;
      for(int i = 0; i < mDim; i++) {
         mW[i].resize(mDim, 0);
         for(int j = 0; j < mDim; j++) {
            int diff = std::min(abs(i-j), mDim - abs(i-j));
            float weight = pow(mHourlyCorr, (float) diff*hoursBetween);
            mW[i][j] = weight * factor;
         }
         // assert(mW[i][i] ==  1);
      }
   }
   return mW;
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

KalmanParameters KalmanFilter::initialize() const {
   std::vector<float> x;
   std::vector<float> k;
   vec2 P;
   P.resize(mDim);
   int hoursBetween = 24 / mDim;
   for(int i = 0; i < P.size(); i++) {
      P[i].resize(mDim, 0);
      for(int j = 0; j < P.size(); j++) {
         int diff = std::min(abs(i-j), mDim - abs(i-j));
         float weight = pow(mHourlyCorr, (float) diff*hoursBetween);
         P[i][j] = weight * mPinit;
      }
   }
   x.resize(mDim, 0);
   k.resize(mDim, 0);
   return KalmanParameters(x, k, P);
}

KalmanParameters::KalmanParameters(std::vector<float> x, std::vector<float> k, vec2 P) {
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
   ss << Util::formatDescription("   elevGradient=-0.0065", "What elevation gradient should be used when interpolating the forecasts to the observation points? In degrees C/m, negative means a cooling with height.") << std::endl;
   return ss.str();
}

KalmanParameters KalmanFilter::update(float iBias, int iTimeStep, const KalmanParameters& iParameters) {
   std::vector<float> x = iParameters.x;
   std::vector<float> k = iParameters.k;
   vec2 P = iParameters.P;

   const vec2& W = getW();

   // Compute Pt|t-1. This increases the covariance.
   for(int i = 0; i < mDim; i++) {
      for(int j = 0; j < mDim; j++) {
         P[i][j] += W[i][j];
      }
   }
   if(Util::isValid(iBias)) {
      float y = iBias;
      int dindex = getParameterIndex(iTimeStep);
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
   }
   else {
      // Missing obs or forecast. P has already been increased, do not update x.
      // TODO: Does the kalman gain need to be updated?
   }

   // Store the new parameters
   KalmanParameters parNew = iParameters;
   parNew.x = x;
   parNew.k = k;
   parNew.P = P;

   return parNew;
}
