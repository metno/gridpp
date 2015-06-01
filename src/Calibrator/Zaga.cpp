#include "Zaga.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile.h"
#include "../Parameters.h"
const float CalibratorZaga::mMaxEnsMean = 100;
CalibratorZaga::CalibratorZaga(const ParameterFile* iParameterFile, Variable::Type iMainPredictor):
      Calibrator(),
      mParameterFile(iParameterFile),
      mMainPredictor(iMainPredictor),
      mFracThreshold(0.5) {
}

bool CalibratorZaga::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   int numInvalidRaw = 0;
   int numInvalidCal = 0;

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Parameters parameters = mParameterFile->getParameters(t);
      Field& precip = *iFile.getField(Variable::Precip, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            std::vector<float> precipRaw = precip(i,j);

            // Compute model variables
            float ensMean = 0;
            float ensFrac = 0;
            int counter = 0;
            bool isValid = true;
            // Check if current ensemble for this gridpoint/time has any missing
            // values. If so, don't calibrate the ensemble.
            for(int e = 0; e < nEns; e++) {
               float value = precipRaw[e];
               if(!Util::isValid(value)) {
                  isValid = false;
                  break;
               }
               ensMean += value;
               ensFrac += (value <= mFracThreshold);
               counter++;
            }
            if(isValid) {
               ensMean = ensMean / counter;
               ensFrac = ensFrac / counter;

               // Limit the input to the calibration model to prevent it
               // from creating very extreme values.
               if(ensMean > mMaxEnsMean) {
                  ensMean = mMaxEnsMean;
               }

               // Calibrate
               std::vector<std::pair<float,int> > pairs(nEns);
               std::vector<float> valuesCal(nEns);
               for(int e = 0; e < nEns; e++) {
                  float quantile = ((float) e+0.5)/nEns;
                  float valueCal   = getInvCdf(quantile, ensMean, ensFrac, parameters);
                  precip(i,j,e) = valueCal;
                  if(!Util::isValid(valueCal))
                     isValid = false;
               }
               if(isValid) {
                  std::vector<float> precipCal = precip(i,j);
                  Calibrator::shuffle(precipRaw, precipCal);
                  for(int e = 0; e < nEns; e++) {
                     precip(i,j,e) = precipCal[e];
                  }
               }
               else {
                  numInvalidCal++;
                  // Calibrator produced some invalid members. Revert to the raw values.
                  for(int e = 0; e < nEns; e++) {
                     precip(i,j,e) = precipRaw[e];
                  }
               }
            }
            else {
               numInvalidRaw++;
               // One or more members are missing, don't calibrate
               for(int e = 0; e < nEns; e++) {
                  precip(i,j,e) = precipRaw[e];
               }
            }
         }
      }
   }
   if(numInvalidRaw > 0) {
      std::stringstream ss;
      ss << "File '" << iFile.getFilename() << "' has " << numInvalidRaw
         << " missing ensembles, out of " << nTime * nLat * nLon << ".";
      Util::warning(ss.str());
   }
   if(numInvalidCal > 0) {
      std::stringstream ss;
      ss << "Calibrator produced " << numInvalidCal
         << " invalid ensembles, out of " << nTime * nLat * nLon << ".";
      Util::warning(ss.str());
   }
   return true;
}

float CalibratorZaga::getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, Parameters& iParameters) {
   if(iQuantile < 0 || iQuantile >= 1) {
      Util::warning("Quantile must be in the interval [0,1)");
      return Util::MV;
   }
   if(!Util::isValid(iEnsMean) || !Util::isValid(iEnsFrac))
      return Util::MV;

   if(iEnsMean < 0 || iEnsFrac < 0 || iEnsFrac > 1)
      return Util::MV;

   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i]))
         return Util::MV;
   }

   if(iQuantile == 0)
      return 0;

   // Check if we are in the discrete mass
   float P0 = getP0(iEnsMean, iEnsFrac, iParameters);
   if(!Util::isValid(P0))
      return Util::MV;
   if(iQuantile < P0)
      return 0;

   float mua = iParameters[0];
   float mub = iParameters[1];
   float sa  = iParameters[2];
   float sb  = iParameters[3];

   float quantileCont = (iQuantile-P0)/(1-P0);
   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = exp(mua + mub * pow(iEnsMean, 1.0/3));
   float sigma = exp(sa + sb * iEnsMean);

   if(mu <= 0 || sigma <= 0)
      return Util::MV;
   if(!Util::isValid(mu) || !Util::isValid(sigma))
      return Util::MV;

   // Parameters in boost and wikipedia
   float shape = 1/(sigma*sigma); // k
   float scale = sigma*sigma*mu;  // theta
   if(!Util::isValid(scale) || !Util::isValid(shape))
      return Util::MV;

   // std::cout << mu << " " << sigma << " " << P0 << " " << shape << " " << scale << std::endl;
   boost::math::gamma_distribution<> dist(shape, scale);
   float value = boost::math::quantile(dist, quantileCont);
   if(!Util::isValid(value))
      return Util::MV;
   return value;
}
float CalibratorZaga::getP0(float iEnsMean, float iEnsFrac, Parameters& iParameters) {
   if(!Util::isValid(iEnsMean) || !Util::isValid(iEnsFrac) || iEnsMean < 0 || iEnsFrac < 0 || iEnsFrac > 1)
      return Util::MV;
   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i]))
         return Util::MV;
   }
   float a = iParameters[4];
   float b = iParameters[5];
   float c = iParameters[6];
   float d = iParameters[7];
   float logit = a + b * iEnsMean + c * iEnsFrac + d * pow(iEnsMean, 1.0/3);
   float P0 = Util::invLogit(logit);
   return P0;
}
void CalibratorZaga::setFracThreshold(float iFraction) {
   if(!Util::isValid(iFraction) || iFraction < 0) {
      std::stringstream ss;
      ss << "CalibratorZaga: fraction threshold (" << iFraction << ") must be 0 or greater.";
      Util::error(ss.str());
   }
   mFracThreshold = iFraction;
}

std::string CalibratorZaga::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c zaga", "Calibrates an ensemble using a zero-adjusted gamma distribution, suitable for parameters like precip and precip. The distribution has three parameters:") << std::endl;
   ss << Util::formatDescription("", "* mean  = exp(a + b * ensmean^(1/3)") << std::endl;
   ss << Util::formatDescription("", "* sigma = exp(c + d * ensmean") << std::endl;
   ss << Util::formatDescription("", "* logit(p0) = e + f * ensmean + g * ensfrac + h * ensmean^(1/3)") << std::endl;
   ss << Util::formatDescription("", "where ensmean is the ensemble mean, and ensfrac is the fraction of members above a certain threshold (use fracThreshold option).") << std::endl;
   ss << Util::formatDescription("   parameters=required", "Read parameters from this text file. The file format is:") << std::endl;
   ss << Util::formatDescription("", "offset0 a b c d e f g h") << std::endl;
   ss << Util::formatDescription("", "...") << std::endl;
   ss << Util::formatDescription("", "offsetN a b c d e f g h") << std::endl;
   ss << Util::formatDescription("", "If the file only has a single line, then the same set of parameters are used for all offsets.") << std::endl;
   ss << Util::formatDescription("   fracThreshold=0.5", "Threshold defining precip/no-precip boundary when computing fraction of members with precip.") << std::endl;
   return ss.str();
}
