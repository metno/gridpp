#include "Zaga.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Parameters.h"
CalibratorZaga::CalibratorZaga(const ParameterFile* iParameterFile, Variable::Type iMainPredictor, const Options& iOptions):
      Calibrator(),
      mParameterFile(iParameterFile),
      mMainPredictor(iMainPredictor),
      mFracThreshold(0.5),
      mOutputPop(false),
      mNeighbourhoodSize(0),
      mPopThreshold(0.5),
      mMaxEnsMean(100),
      m6h(false) {
   iOptions.getValue("fracThreshold", mFracThreshold);
   iOptions.getValue("outputPop", mOutputPop);
   iOptions.getValue("neighbourhoodSize", mNeighbourhoodSize);
   iOptions.getValue("popThreshold", mPopThreshold);
   iOptions.getValue("maxEnsMean", mMaxEnsMean);
   iOptions.getValue("6h", m6h);
   if(!Util::isValid(mMaxEnsMean) || mMaxEnsMean <= 0) {
      Util::error("CalibratorZaga maxEnsMean must be > 0");
   }
   if(mNeighbourhoodSize < 0) {
      std::stringstream ss;
      ss << "CalibratorZaga neighbourhoodSize (" << mNeighbourhoodSize << ") must be >= 0";
      Util::error(ss.str());
   }
}

bool CalibratorZaga::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   Variable::Type popVariable = Variable::Pop;
   int startTime = 0;
   int timeWindow = 1;
   if(m6h) {
      popVariable = Variable::Pop6h;
      startTime = 5;
      timeWindow = 6;
   }

   std::vector<FieldPtr> precips;
   for(int t = 0; t < nTime; t++) {
      precips.push_back(iFile.getField(Variable::Precip, t));
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      int numInvalidRaw = 0;
      int numInvalidCal = 0;

      Parameters parameters = mParameterFile->getParameters(t);

      // Load the POP output field, if needed
      FieldPtr pop;
      if(mOutputPop) {
         pop = iFile.getField(popVariable, t);
      }

      #pragma omp parallel for reduction(+:numInvalidRaw, numInvalidCal)
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            // for Pop6h, the first few hours are undefined, since we cannot do a 6h accumulation
            if(mOutputPop && t < startTime) {
               for(int e = 0; e < nEns; e++) {
                  (*pop)(i,j,e) = Util::MV;
               }
            }
            else {
               // Compute model variables
               float ensMean = 0;
               float ensFrac = 0;
               int counter = 0;
               bool isValid = true;
               for(int e = 0; e < nEns; e++) {
                  // Create a neighbourhood ensemble
                  for(int ii = std::max(0, i-mNeighbourhoodSize); ii <= std::min(nLat-1, i+mNeighbourhoodSize); ii++) {
                     for(int jj = std::max(0, j-mNeighbourhoodSize); jj <= std::min(nLon-1, j+mNeighbourhoodSize); jj++) {

                        float total = 0;
                        // Accumulate over the time window, if necessary
                        for(int tt = t-timeWindow + 1; tt <= t; tt++) {
                           Field& precip = *precips[tt];

                           float value = precip(ii,jj,e);
                           if(Util::isValid(total) && Util::isValid(value)) {
                              total += value;
                           }
                           else {
                              total = Util::MV;
                           }
                        }
                        if(Util::isValid(total) && Util::isValid(ensMean)) {
                           ensMean += total;
                           ensFrac += (total <= mFracThreshold);
                           counter++;
                        }
                        else {
                           ensMean = Util::MV;
                           ensFrac = Util::MV;
                        }
                     }
                  }
               }
               Field& precip = *precips[t];
               const std::vector<float>& precipRaw = precip(i,j);

               // Only calibrate the ensemble if all members are available. Otherwise
               // use the raw members.
               if(Util::isValid(ensMean) && counter > 0) {
                  ensMean = ensMean / counter;
                  ensFrac = ensFrac / counter;

                  if(i == 0 && j == 0) {
                     // std::cout << ensMean << " " << ensFrac << " " << counter << std::endl;
                  }

                  // Limit the input to the calibration model to prevent it
                  // from creating very extreme values.
                  if(Util::isValid(mMaxEnsMean) && ensMean > mMaxEnsMean) {
                     ensMean = mMaxEnsMean;
                  }

                  // Compute POP
                  if(mOutputPop) {
                     for(int e = 0; e < nEns; e++) {
                        float cdf = getCdf(mPopThreshold, ensMean, ensFrac, parameters);
                        // std::cout << ensMean << " " << ensFrac << " " << cdf << std::endl;
                        if(Util::isValid(cdf))
                           (*pop)(i,j,e) = 1 - cdf;
                        else
                           (*pop)(i,j,e) = Util::MV;
                     }
                  }
                  // Compute precip ensemble
                  else {
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
                        // Calibrator produced some invalid members. Revert to the raw values.
                        for(int e = 0; e < nEns; e++) {
                           precip(i,j,e) = precipRaw[e];
                        }
                     }
                  }
               }
               else {
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
         ss << "File '" << iFile.getFilename() << "' missing " << numInvalidRaw
            << "/" << nLat * nLon << " ensembles for timestep " << t << ".";
         Util::warning(ss.str());
      }
      if(numInvalidCal > 0) {
         std::stringstream ss;
         ss << "Calibrator produces '" << numInvalidRaw
            << "/" << nLat * nLon << " invalid ensembles for timestep " << t << ".";
         Util::warning(ss.str());
      }
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
float CalibratorZaga::getCdf(float iThreshold, float iEnsMean, float iEnsFrac, Parameters& iParameters) {
   if(!Util::isValid(iThreshold) || !Util::isValid(iEnsMean) || !Util::isValid(iEnsFrac))
      return Util::MV;

   if(iEnsMean < 0 || iEnsFrac < 0 || iEnsFrac > 1)
      return Util::MV;

   if(iThreshold < 0)
      return 0;

   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i]))
         return Util::MV;
   }

   // Check if we are in the discrete mass
   float P0 = getP0(iEnsMean, iEnsFrac, iParameters);
   if(!Util::isValid(P0))
      return Util::MV;
   if(iThreshold == 0)
      return P0;

   float mua = iParameters[0];
   float mub = iParameters[1];
   float sa  = iParameters[2];
   float sb  = iParameters[3];

   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = exp(mua + mub * pow(iEnsMean, 1.0/3));
   float sigma = exp(sa + sb * iEnsMean);

   // std::cout << mua << " " << mub << " " << sa << " " << sb << " " << mu << " " << sigma << std::endl;

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
   float contCdf = boost::math::cdf(dist, iThreshold) ;
   float cdf = P0 + (1 - P0)*contCdf;
   // std::cout << cdf << " " << P0 << " " << iThreshold << " " << contCdf << " " << shape << " " << scale << std::endl;
   if(!Util::isValid(cdf))
      return Util::MV;
   assert(cdf <= 1);
   assert(cdf >= 0);
   return cdf;
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
   ss << Util::formatDescription("   neighbourhoodSize=0", "Increase the ensemble by taking all gridpoints within a neighbourhood. A value of 0 means no neighbourhood is used.") << std::endl;
   ss << Util::formatDescription("   outputPop=0", "Should probability of precip be written to the POP field?") << std::endl;
   ss << Util::formatDescription("   popThreshold=0.5", "If POP is written, what threshold should be used?") << std::endl;
   ss << Util::formatDescription("   maxEnsMean=100", "Upper limit of what the ensemble mean is allowed to be when passed into the distribution. This effectively prevents the distribution to yield very high values.") << std::endl;
   ss << Util::formatDescription("   6h=0", "If POP is produced, should it be based on the precip in the last 6 hours? If so, the Pop6h variable is written.") << std::endl;
   return ss.str();
}
