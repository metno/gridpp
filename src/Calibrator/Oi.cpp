#include "Oi.h"
#include "../Util.h"
#include "../Parameters.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <math.h>
#include <armadillo>
#include "/home/thomasn/local/include/armadillo"

CalibratorOi::CalibratorOi(Variable iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mVLength(100),
      mHLength(30000),
      mMu(0.9),
      mMinObs(3),
      mSort(true),
      mMinRho(0.0013),
      mEpsilon(0.5),
      mUseRho(true),
      mObsOnly(false),
      mElevGradient(-0.0065),
      mMethod(Util::MV),
      mBiasVariable(""),
      mSigma(1),
      mC(1.03),
      mSaveDiff(false),
      mDeltaVariable(""),
      // Add mDeltaVariable
      mX(Util::MV),
      mY(Util::MV),
      mMaxLocations(20),
      mNumVariable(""),
      mUseMeanBias(false),
      mMaxElevDiff(200),
      // Default model error variance
      mMinValidEns(5),
      mTest(Util::MV),
      mNewDeltaVar(1),
      mDiagnose(false),
      mWMin(0.5),
      mMaxBytes(6.0 * 1024 * 1024 * 1024),
      mGamma(0.25) {
   iOptions.getValue("biasVariable", mBiasVariable);
   iOptions.getValue("d", mHLength);
   iOptions.getValue("h", mVLength);
   iOptions.getValue("maxLocations", mMaxLocations);
   iOptions.getValue("sort", mSort);
   iOptions.getValue("sigma", mSigma);
   iOptions.getValue("deltaVariable", mDeltaVariable);
   iOptions.getValue("gamma", mGamma);
   iOptions.getValue("obsOnly", mObsOnly);
   iOptions.getValue("mu", mMu);
   iOptions.getValue("minObs", mMinObs);
   iOptions.getValue("x", mX);
   iOptions.getValue("y", mY);
   iOptions.getValue("minRho", mMinRho);
   iOptions.getValue("useRho", mUseRho);
   iOptions.getValue("saveDiff", mSaveDiff);
   iOptions.getValue("maxBytes", mMaxBytes);
   iOptions.getValue("method", mMethod);
   iOptions.getValue("minEns", mMinValidEns);
   iOptions.getValue("numVariable", mNumVariable);
   iOptions.getValue("elevGradient", mElevGradient);
   iOptions.getValue("useMeanBias", mUseMeanBias);
   iOptions.getValue("wmin", mWMin);
   iOptions.getValue("epsilon", mEpsilon);
   iOptions.getValue("test", mTest);
   iOptions.getValue("c", mC);
   iOptions.getValue("diagnose", mDiagnose);
   iOptions.getValue("newDeltaVar", mNewDeltaVar);
   iOptions.getValue("maxElevDiff", mMaxElevDiff);  // Don't use obs that are further than this from their nearest neighbour
   iOptions.check();

   // Gamma: The error covariance of the bias is this fraction of the background error
}

// Set up convenient functions for debugging in gdb
template<class Matrix>
void print_matrix(Matrix matrix) {
       matrix.print(std::cout);
}
typedef arma::mat mattype;
typedef arma::vec vectype;
typedef arma::cx_mat cxtype;
template void print_matrix<mattype>(mattype matrix);
template void print_matrix<cxtype>(cxtype matrix);

bool CalibratorOi::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();
   vec2 lafs = iFile.getLandFractions();

   // Check if this method can be applied
   bool hasValidGridpoint = false;
   for(int y = 0; y < nY; y++) {
      for(int x = 0; x < nX; x++) {
         if(Util::isValid(lats[y][x]) && Util::isValid(lons[y][x]) && Util::isValid(elevs[y][x])) {
            hasValidGridpoint = true;
         }
      }
   }
   if(!hasValidGridpoint) {
      Util::warning("There are no gridpoints with valid lat/lon/elev values. Skipping oi...");
      return false;
   }

   std::vector<Location> obsLocations = iParameterFile->getLocations();
   if(iParameterFile->getNumParameters() != 2) {
      std::stringstream ss;
      ss << "Parameter file has " << iParameterFile->getNumParameters() << " parameters, not 2";
      Util::error(ss.str());
   }

   // Retrieve parameters
   std::vector<float> ci(obsLocations.size());
   std::vector<float> obs(obsLocations.size());
   std::vector<float> obselevs(obsLocations.size());
   for(int i = 0; i < obsLocations.size(); i++) {
      Parameters parameters = iParameterFile->getParameters(0, obsLocations[i]);
      obs[i] = parameters[0];
      obselevs[i] = obsLocations[i].elev();
      ci[i] = parameters[1];
   }
   float sigma = mSigma;
   int S = obsLocations.size();
   // Find the spacing between each grid
   float gridSize = Util::getDistance(lats[0][0], lons[0][0], lats[1][0], lons[1][0]);
   std::stringstream ss;
   ss << "Grid size: " << gridSize << " m";
   Util::info(ss.str());

   // Loop over each observation, find the nearest gridpoint and place the obs into all gridpoints
   // in the vicinity of the nearest neighbour. This is only meant to be an approximation, but saves
   // considerable time instead of doing a loop over each grid point and each observation.

   // Store the indicies (into the obsLocations array) that a gridpoint has available
   std::vector<std::vector<std::vector<int> > > obsIndices; // Y, X, obs indices
   std::vector<float> obsY(S);
   std::vector<float> obsX(S);
   std::vector<float> obsLaf(S);
   obsIndices.resize(nY);
   for(int y = 0; y < nY; y++) {
      obsIndices[y].resize(nX);
   }

   // Spread each observation out to this many gridpoints from the nearest neighbour
   int radius = 3.64 * mHLength / gridSize;

   // When large radiuses are used, the process becomes memory-intensive. Try to fail here
   // if we expect to use more memory than desired. The true memory is roughly
   // 1 GB + expectedBytes * F
   int bytesPerValue = 4;
   float expectedBytes = float(radius * radius) * 4 * bytesPerValue * S;
   std::cout << "Expected MB: " << 1000 + expectedBytes / 1024 / 1024 << std::endl;
   if(Util::isValid(mMaxBytes) && expectedBytes > mMaxBytes) {
      std::stringstream ss;
      ss << "Expected size (" << expectedBytes / 1024 / 1024 << " GB) is greater than "
         << float(mMaxBytes) / 1024 / 1024 << " GB";
      Util::error(ss.str());
   }

   double time_s = Util::clock();
   KDTree searchTree(iFile.getLats(), iFile.getLons());
   int count = 0;
   for(int i = 0; i < S; i++) {
      if(i % 1000 == 0) {
         std::stringstream ss;
         ss << i;
         Util::info(ss.str());
      }
      if(Util::isValid(obs[i])) {
         int Y, X;
         searchTree.getNearestNeighbour(obsLocations[i].lat(), obsLocations[i].lon(), Y, X);
         // Check if the elevation of the station against the reference grid
         bool hasValidElev = Util::isValid(obsLocations[i].elev());
         if(Util::isValid(mMaxElevDiff) && hasValidElev) {
            float elevDiff = abs(obsLocations[i].elev() - elevs[Y][X]);
            hasValidElev = elevDiff < mMaxElevDiff;
         }
         if(hasValidElev) {
            for(int y = std::max(0, Y - radius); y < std::min(nY, Y + radius); y++) {
               for(int x = std::max(0, X - radius); x < std::min(nX, X + radius); x++) {
                  if(mSort || obsIndices[y][x].size() < mMaxLocations) {
                     obsIndices[y][x].push_back(i);
                     count ++;
                  }
               }
            }
         }
         else {
            std::stringstream ss;
            ss << "Removing station because elevation (" << obsLocations[i].elev() << " m) is too far from grid (" << elevs[Y][X] << " m)";
            Util::warning(ss.str());
         }
         obsY[i] = Y;
         obsX[i] = X;
         obsLaf[i] = lafs[Y][X];
      }
   }
   double time_e = Util::clock();
   std::cout << "Assigning locations " << time_e - time_s << std::endl;


   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      FieldPtr output = iFile.getEmptyField();
      FieldPtr bias;
      FieldPtr newbias;
      FieldPtr delta;
      FieldPtr newdelta;
      bool useBias = mBiasVariable != "";
      if(useBias) {
         bias = iFile.getField(mBiasVariable, t);

         // Scale the bias
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               if(Util::isValid((*bias)(y, x, 0)))
                  (*bias)(y, x, 0) *= mMu;
               else
                  (*bias)(y, x, 0) = 0;
            }
         }
         newbias = iFile.getEmptyField(0);
      }
      bool useDelta = mDeltaVariable != "";
      if(useDelta) {
         delta = iFile.getField(mDeltaVariable, t);
         // Initialize if missing
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               if(!Util::isValid((*delta)(y, x, 0)))
                  (*delta)(y, x, 0) = 1;
            }
         }
         newdelta = iFile.getEmptyField(0);
      }

      FieldPtr num;
      if(mNumVariable != "")
         num = iFile.getField(mNumVariable, t);

      // Compute Y
      vec2 Yglobal(S);
      std::vector<float> Yhatglobal(S);
      for(int i = 0; i < S; i++) {
         Yglobal[i].resize(nEns, 0);
         float elevCorr = 0;
         if(Util::isValid(mElevGradient)) {
            float nnElev = elevs[obsY[i]][obsX[i]];
            assert(Util::isValid(nnElev));
            assert(Util::isValid(obsLocations[i].elev()));
            float elevDiff = obsLocations[i].elev() - nnElev;
            elevCorr = mElevGradient * elevDiff;
         }
         float total = 0;
         int count = 0;
         for(int e = 0; e < nEns; e++) {
            float value = (*field)(obsY[i], obsX[i], e);
            if(Util::isValid(value)) {
               value += elevCorr;
               Yglobal[i][e] = value;
               total += value;
               count++;
            }
         }
         float mean = total / count;
         for(int e = 0; e < nEns; e++) {
            float value = Yglobal[i][e];
            if(Util::isValid(value)) {
               Yglobal[i][e] -= mean;
            }
         }
         Yhatglobal[i] = mean;
         if(useBias) {
            float currBias = (*bias)(obsY[i], obsX[i], 0);
            Yhatglobal[i] -= currBias;
         }
      }

      #pragma omp parallel for
      for(int x = 0; x < nX; x++) {
         for(int y = 0; y < nY; y++) {
            float lat = lats[y][x];
            float lon = lons[y][x];
            float elev = elevs[y][x];
            float laf = lafs[y][x];
            std::vector<int> useLocations0 = obsIndices[y][x];

            //
            // Create list of locations for this gridpoint
            //
            std::vector<int> useLocations;
            useLocations.reserve(useLocations0.size());
            std::vector<std::pair<float,int> > rhos0;
            rhos0.reserve(useLocations0.size());
            for(int i = 0; i < useLocations0.size(); i++) {
               int index = useLocations0[i];
               float hdist = Util::getDistance(obsLocations[index].lat(), obsLocations[index].lon(), lat, lon, true);
               float vdist = obsLocations[index].elev() - elev;
               float lafdist = 1;
               if(Util::isValid(obsLaf[i]) && Util::isValid(laf))
                  lafdist = obsLaf[index] - laf;
               float rho = calcRho(hdist, vdist, lafdist);
               int X = obsX[index];
               int Y = obsY[index];
               // Only include observations that are within the domain
               if(X > 0 && X < lats[0].size()-1 && Y > 0 && Y < lats.size()-1) {
                  if(rho > mMinRho) {
                     rhos0.push_back(std::pair<float,int>(rho, i));
                  }
               }
            }
            arma::vec rhos;
            if(rhos0.size() > mMaxLocations) {
               // If sorting is enabled and we have too many locations, then only keep the best ones based on rho.
               // Otherwise, just use the last locations added
               rhos = arma::vec(mMaxLocations);
               if(mSort) {
                  std::sort(rhos0.begin(), rhos0.end(), Util::sort_pair_first<float,int>());
               }
               for(int i = 0; i < mMaxLocations; i++) {
                  // The best values start at the end of the array
                  int index = rhos0[rhos0.size() - 1 - i].second;
                  useLocations.push_back(useLocations0[index]);
                  rhos(i) = rhos0[rhos0.size() - 1 - i].first;
               }
            }
            else {
               rhos = arma::vec(rhos0.size());
               for(int i = 0; i < rhos0.size(); i++) {
                  int index = rhos0[i].second;
                  useLocations.push_back(useLocations0[index]);
                  rhos(i) = rhos0[i].first;
               }
            }

            int nObs = useLocations.size();

            if(nObs < mMinObs) {
               // If we have too few observations though, then use the background
               for(int e = 0; e < nEns; e++) {
                  if(mSaveDiff)
                     (*output)(y, x, e) = Util::MV;
                  else
                     (*output)(y, x, e) = (*field)(y, x, e);
               }
               continue;
            }

            if(mObsOnly) {
               // Here we don't run the OI algorithm but instead just use the median
               // of all observations
               std::vector<float> currObs(nObs, Util::MV);
               for(int i = 0; i < nObs; i++) {
                  int index = useLocations[i];
                  currObs[i] = obs[index];
               }
               float value = Util::calculateStat(currObs, Util::StatTypeQuantile, 0.5);
               for(int e = 0; e < nEns; e++) {
                  if(mSaveDiff)
                     (*output)(y, x, e) = value - (*field)(y, x, e);
                  else
                     (*output)(y, x, e) = value;
               }
               continue;
            }

            int nValidEns = 0;
            std::vector<int> validEns;
            for(int e = 0; e < nEns; e++) {
               float value = (*field)(y, x, e);
               if(Util::isValid(value)) {
                  validEns.push_back(e);
                  nValidEns++;
               }
            }
            if(nValidEns < mMinValidEns) {
               for(int e = 0; e < nEns; e++) {
                  (*output)(y, x, e) = (*field)(y, x, e);
               }
               continue;
            }

            vectype currObs(nObs);
            vectype currElev(nObs);
            vectype currLaf(nObs);
            for(int i = 0; i < useLocations.size(); i++) {
               int index = useLocations[i];
               float curr = obs[index];
               float elev = obselevs[index];
               float laf = obsLaf[index];
               currObs[i] = curr;
               currElev[i] = elev;
               currLaf[i] = laf;
            }

            // Compute Rinv
            mattype Rinv(nObs, nObs, arma::fill::zeros);
            for(int i = 0; i < nObs; i++) {
               int index = useLocations[i];
               float r = sigma * sigma * ci[index];
               if(mUseRho) {
                  float rho = rhos[i];
                  Rinv(i, i) = 1 / r * rho;
               }
               else {
                  Rinv(i, i) = 1 / r;
               }
            }

            // Compute Y (model at obs-locations)
            mattype Y(nObs, nValidEns);
            vectype Yhat(nObs);

            for(int i = 0; i < nObs; i++) {
               // Use the nearest neighbour for this location
               int index = useLocations[i];
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  Y(i, e) = Yglobal[index][ei];
               }
               Yhat(i) = Yhatglobal[index];

            }

            // Compute C matrix
            // k x S * S x S
            mattype C(nValidEns, nObs);
            C = Y.t() * Rinv;

            mattype Pinv(nValidEns, nValidEns);
            float currDelta = 1;
            if(useDelta) {
               currDelta = (*delta)(y, x, 0);
            }
            float diag = 1 / currDelta * (nValidEns - 1);
            if(useBias)
               diag = 1 / currDelta / (1 + mGamma) * (nValidEns - 1);

            Pinv = C * Y + diag * arma::eye<mattype>(nValidEns, nValidEns);
            float cond = arma::rcond(Pinv);
            if(cond <= 0) {
               std::stringstream ss;
               ss << "Condition number of " << cond << ". Using raw values";
               Util::warning(ss.str());
               for(int e = 0; e < nEns; e++) {
                  (*output)(y, x, e) = (*field)(y, x, e); // Util::MV;
               }
               continue;
            }

            mattype P = arma::inv(Pinv);
            cxtype Wcx(nValidEns, nValidEns);
            bool status = arma::sqrtmat(Wcx, (nValidEns - 1) * P);
            if(!status) {
               std::cout << "Near singular matrix for sqrtmat:" << std::endl;
               std::cout << "Lat: " << lat << std::endl;
               std::cout << "Lon: " << lon << std::endl;
               std::cout << "Elev: " << elev << std::endl;
               std::cout << "Laf: " << laf << std::endl;
               std::cout << "Pinv" << std::endl;
               print_matrix<mattype>(Pinv);
               std::cout << "P" << std::endl;
               print_matrix<mattype>(P);
               std::cout << "Y:" << std::endl;
               print_matrix<mattype>(Y);
               std::cout << "currObs:" << std::endl;
               print_matrix<mattype>(currObs);
               std::cout << "Yhat" << std::endl;
               print_matrix<mattype>(Yhat);
            }

            mattype W = arma::real(Wcx);
            if(W.n_rows == 0) {
               std::stringstream ss;
               ss << "Could not find the real part of W. Using raw values.";
               Util::warning(ss.str());
               for(int e = 0; e < nEns; e++) {
                  (*output)(y, x, e) = (*field)(y, x, e);
               }
               continue;
            }

            // Compute PC
            mattype PC(nValidEns, nObs);
            PC = P * C;

            // Compute w
            vectype w(nValidEns);
            if(mDiagnose)
               w = PC * (arma::ones<vectype>(nObs));
            else
               w = PC * (currObs - Yhat);

            // Add w to W. TODO: Is this done correctly?
            for(int e = 0; e < nValidEns; e++) {
               for(int e2 = 0; e2 < nValidEns; e2 ++) {
                  W(e, e2) = W(e, e2) + w(e) ;
               }
            }

            // Compute X (perturbations about model mean)
            vectype X(nValidEns);
            float total = 0;
            int count = 0;
            for(int e = 0; e < nValidEns; e++) {
               int ei = validEns[e];
               float value = (*field)(y, x, ei);
               if(Util::isValid(value)) {
                  X(e) = value;
                  total += value;
                  count++;
               }
               else {
                  std::cout << "Invalid value " << y << " " << x << " " << e << std::endl;
               }
            }
            float mean = total / count;
            for(int e = 0; e < nValidEns; e++) {
               X(e) -= mean;
            }

            // Write debugging information
            if(x == mX && y == mY) {
               std::cout << "Lat: " << lat << std::endl;
               std::cout << "Lon: " << lon << std::endl;
               std::cout << "Elev: " << elev << std::endl;
               std::cout << "Laf: " << laf << std::endl;
               std::cout << "Num obs: " << nObs << std::endl;
               std::cout << "Num ens: " << nValidEns << std::endl;
               std::cout << "rhos" << std::endl;
               print_matrix<mattype>(rhos);
               std::cout << "P" << std::endl;
               print_matrix<mattype>(P);
               std::cout << "C" << std::endl;
               print_matrix<mattype>(C);
               std::cout << "PC" << std::endl;
               print_matrix<mattype>(PC);
               std::cout << "W" << std::endl;
               print_matrix<mattype>(W);
               std::cout << "w" << std::endl;
               print_matrix<mattype>(w);
               std::cout << "Y:" << std::endl;
               print_matrix<mattype>(Y);
               std::cout << "Yhat" << std::endl;
               print_matrix<mattype>(Yhat);
               std::cout << "currObs" << std::endl;
               print_matrix<mattype>(currObs);
               std::cout << "currObs - Yhat" << std::endl;
               print_matrix<mattype>(currObs - Yhat);
               std::cout << "X" << std::endl;
               print_matrix<mattype>(X);
               std::cout << "elevs" << std::endl;
               print_matrix<mattype>(currElev);
               std::cout << "lafs" << std::endl;
               print_matrix<mattype>(currLaf);
               std::cout << "Analaysis increment:" << std::endl;
               print_matrix<mattype>(X.t() * W);
               std::cout << "My: " << arma::mean(arma::dot(currObs - Yhat, rhos) / nObs) << std::endl;
            }

            // Compute analysis
            for(int e = 0; e < nValidEns; e++) {
               int ei = validEns[e];
               float total = 0;
               for(int k = 0; k < nValidEns; k++) {
                  total += X(k) * W(k, e);
               }

               float diff = total;
               if(mUseMeanBias)
                  diff = arma::mean(currObs - Yhat);

               if(mSaveDiff)
                  (*output)(y, x, ei) = diff;
               else {
                  float raw = (*field)(y, x, ei);
                  if(useBias) {
                     raw -= (*bias)(y, x, 0);
                  }
                  (*output)(y, x, ei) = raw + diff;
               }

               if(mNumVariable != "")
                  (*num)(y, x, ei) = nObs;

            }
            // Update bias
            if(useBias) {
               float biasTotal = 0;
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  biasTotal += (*field)(y, x, ei) * w(e);
               }
               (*newbias)(y, x, 0) = (*bias)(y, x, 0) - mGamma / (1 + mGamma) * biasTotal;
            }


            // Update delta
            /*
            float deltaVar = mC - 1;
            float trace = arma::trace(Y * Y.t());
            float numerator = mSigma * mSigma / mEpsilon / mEpsilon;
            float denomenator = 1.0 / nObs / (nValidEns - 1) * trace;
            float currDeltaEvidence = numerator / denomenator;
            float weightOld = deltaVar;
            float weightNew = mNewDeltaVar;
            (*newdelta)(y, x, 0) = ((*delta)(y, x, 0) * weightNew + currDeltaEvidence * weightOld) / (weightOld + weightNew);
            */
         }
      }
      std::cout << "Adding" << std::endl;
      iFile.addField(output, mVariable, t);
      if(mNumVariable != "")
         iFile.addField(num, Variable(mNumVariable), t);
      if(useBias) {
         std::cout << "Adding bias field" << std::endl;
         iFile.addField(newbias, Variable(mBiasVariable), t);
      }
      if(useDelta) {
         // Update delta
         float deltaVar = mC - 1;
         float trace = 0;
         float numValidS = 0;
         for(int s = 0; s < S; s++) {
            // Compute value in coord s,s
            float value = 0;
            int count = 0;
            for(int e = 0; e < nEns; e++) {
               if(Util::isValid(Yglobal[s][e])) {
                  value += Yglobal[s][e] * Yglobal[s][e];
                  count++;
               }
            }
            if(count > 1) {
               value = value / (count-1);
               trace += value;
               numValidS++;
            }
         }
         float numerator = mSigma * mSigma / mEpsilon / mEpsilon;
         float denomenator = 1.0 / numValidS * trace;
         float currDeltaEvidence = numerator / denomenator;
         float weightOld = deltaVar;
         float weightNew = mNewDeltaVar;
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               (*newdelta)(y, x, 0) = ((*delta)(y, x, 0) * weightNew + currDeltaEvidence * weightOld) / (weightOld + weightNew);
            }
         }

         std::cout << "Adding delta field" << std::endl;
         iFile.addField(newdelta, Variable(mDeltaVariable), t);
      }
   }
   return true;
}

float CalibratorOi::calcRho(float iHDist, float iVDist, float iLDist) const {
   float h = (iHDist/mHLength);
   float v = (iVDist/mVLength);
   float l = 1 - (1 - mWMin) * std::abs(iLDist);
   float rho = exp(-0.5 * h * h) * exp(-0.5 * v * v) * l;
   return rho;
}

std::string CalibratorOi::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c oi","Spreads bias in space by using kriging. A parameter file is required, which must have one column with the bias.")<< std::endl;
   return ss.str();
}
