#ifndef KALMANFILTER_H
#define KALMANFILTER_H
#include "Variable.h"
class Obs;
class Forecast;
class Parameters;
class File;
class Options;
class ParameterFile;

typedef std::vector<std::vector<float> > vec2;

// Container for KF parameters. Stores:
// - x (best estimate of bias)
// - k (Kalman gain)
// - P (error covariance of best estimate)
class KalmanParameters {
   public:
      KalmanParameters(std::vector<float> x, std::vector<float> k, vec2 P);
      KalmanParameters(const Parameters& iParameters);

      std::vector<float> x;
      std::vector<float> k;
      vec2 P;

      // Serialize parameters into regular gridpp parameters
      Parameters toParameters() const;

      // Dimension of filter
      int dim() const {return x.size();};
};

class KalmanFilter {
   public:
      KalmanFilter(const Variable& iVariable, const Options& iOptions);
      // Run the KF on forecast and observation files for a specified timestep
      // @param iDbIn read KF parameters from this file. If NULL, initialize new ones.
      // @param iDbOut If not NULL, write KF parameters to this file.
      // @param iBiasFile If not NULL, write biases to this file.
      bool writeBiasFile(const File& iFcst,
            const File& iObs,
            int iTimeStep,
            const ParameterFile* iDbIn=NULL,
            ParameterFile* iDbOut=NULL,
            ParameterFile* iBiasFile=NULL);

      // Update the kalman filter parameters based on a measured bias at some time
      KalmanParameters update(float iBias, int iTime, const KalmanParameters& iParameters);

      // Returns KF parameters initialized when no information is available
      KalmanParameters initialize() const;
      static std::string description();
   private:
      vec2 getW() const;
      int getParameterIndex(int iTime) const;

      Variable mVariable;
      float mRatio; // std W / std V
      float mHourlyCorr; // Auto-correlation of bias measurements
      int mDim; // Number of observations to use per 24 hours
      float mElevGradient;
      float mPinit;    // Initial standard error of estimate
      mutable vec2 mW; // Error covariance matrix
      float mV;        // Standard error of bias measurements
};
#endif
