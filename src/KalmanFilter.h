#ifndef KALMANFILTER_H
#define KALMANFILTER_H
#include "Variable.h"
class Obs;
class Forecast;
class Parameters;
class File;
class Options;
class ParameterFile;

class KalmanParameters {
   public:
      // Use existing parameters
      KalmanParameters(const Parameters& iParameters);
      // Initialize new parameters
      KalmanParameters(std::vector<float> x, std::vector<float> k, std::vector<std::vector<float> > P);
      Parameters toParameters() const;
      std::vector<float> x; // Last estimate
      std::vector<float> k; // Kalman gain
      std::vector<std::vector<float> > P; // Covariance matrix
      int dim() const {return x.size();};
};

class KalmanFilter {
   public:
      KalmanFilter(Variable::Type iVariable, const Options& iOptions);
      bool writeBiasFile(const File& iFcst,
            const File& iObs,
            int iCurrDate,
            int iStartTime,
            int iEndTime,
            const ParameterFile* iDbIn=NULL,
            ParameterFile* iDbOut=NULL,
            ParameterFile* iBiasFile=NULL);
      static std::string description();
   private:
      Variable::Type mVariable;
      float mRatio; // std W / std V
      float mHourlyCorr;
      float mV;
      int mDim; // Hours between obs?
      std::vector<std::vector<float> > getW() const;
      int getParameterIndex(int iTime) const;
      KalmanParameters initialize() const;
};
class KalmanDb {
   public:
      KalmanDb(const ParameterFile& iFile);
};
#endif
