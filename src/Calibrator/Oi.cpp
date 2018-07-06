#include "Oi.h"
#include "../Util.h"
#include "../Parameters.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <math.h>
#include <armadillo>

CalibratorOi::CalibratorOi(Variable iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mVLength(100),
      mHLength(30000),
      mMu(0.9),
      mMinObs(0),
      mSort(true),
      mMinRho(0.0013),
      mEpsilon(0.5),
      mUseRho(true),
      mObsOnly(false),
      mElevGradient(-0.0065),
      mMethod(Util::MV),
      mUseBug(false),
      mBiasVariable(""),
      mSigma(1),
      mDelta(-999),
      mC(1.03),
      mSaveDiff(false),
      mDeltaVariable(""),
      mUseEns(true),
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
      mExtrapolate(false),
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
   iOptions.getValue("delta", mDelta);
   iOptions.getValue("deltaVariable", mDeltaVariable);
   iOptions.getValue("gamma", mGamma);
   iOptions.getValue("obsOnly", mObsOnly);
   iOptions.getValue("mu", mMu);
   iOptions.getValue("minObs", mMinObs);
   iOptions.getValue("x", mX);
   iOptions.getValue("y", mY);
   iOptions.getValue("extrapolate", mExtrapolate);
   iOptions.getValue("minRho", mMinRho);
   iOptions.getValue("useRho", mUseRho);
   iOptions.getValue("useBug", mUseBug);
   iOptions.getValue("saveDiff", mSaveDiff);
   iOptions.getValue("maxBytes", mMaxBytes);
   iOptions.getValue("method", mMethod);
   iOptions.getValue("minEns", mMinValidEns);
   iOptions.getValue("numVariable", mNumVariable);
   iOptions.getValue("elevGradient", mElevGradient);
   iOptions.getValue("useMeanBias", mUseMeanBias);
   iOptions.getValue("useEns", mUseEns);
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

template void print_matrix<CalibratorOi::mattype>(CalibratorOi::mattype matrix);
template void print_matrix<CalibratorOi::cxtype>(CalibratorOi::cxtype matrix);

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

   std::vector<Location> gLocations = iParameterFile->getLocations();
   if(iParameterFile->getNumParameters() != 2) {
      std::stringstream ss;
      ss << "Parameter file has " << iParameterFile->getNumParameters() << " parameters, not 2";
      Util::error(ss.str());
   }

   float sigma = mSigma;
   int gS = gLocations.size();
   // Find the spacing between each grid
   float gridSize = Util::getDistance(lats[0][0], lons[0][0], lats[1][0], lons[1][0]);
   std::stringstream ss;
   ss << "Grid size: " << gridSize << " m";
   Util::info(ss.str());

   // Loop over each observation, find the nearest gridpoint and place the obs into all gridpoints
   // in the vicinity of the nearest neighbour. This is only meant to be an approximation, but saves
   // considerable time instead of doing a loop over each grid point and each observation.

   // Store the indicies (into the gLocations array) that a gridpoint has available
   std::vector<std::vector<std::vector<int> > > gLocIndices; // Y, X, obs indices
   std::vector<float> gYi(gS, Util::MV);
   std::vector<float> gXi(gS, Util::MV);
   std::vector<float> gLafs(gS, Util::MV);
   std::vector<float> gElevs(gS, Util::MV);
   std::vector<float> gCi(gS, Util::MV);
   std::vector<float> gObs(gS, Util::MV);
   gLocIndices.resize(nY);
   for(int y = 0; y < nY; y++) {
      gLocIndices[y].resize(nX);
   }

   // Spread each observation out to this many gridpoints from the nearest neighbour
   int radius = 3.64 * mHLength / gridSize;

   // When large radiuses are used, the process becomes memory-intensive. Try to fail here
   // if we expect to use more memory than desired. The true memory is roughly
   // 1 GB + expectedBytes * F
   int bytesPerValue = 4;
   float expectedBytes = float(radius * radius) * 4 * bytesPerValue * gS;
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
   for(int i = 0; i < gS; i++) {
      if(i % 1000 == 0) {
         std::stringstream ss;
         ss << i;
         Util::info(ss.str());
      }
      Parameters parameters = iParameterFile->getParameters(0, gLocations[i]);
      gObs[i] = parameters[0];
      gCi[i] = parameters[1];
      gElevs[i] = gLocations[i].elev();
      if(Util::isValid(gObs[i])) {
         int Y, X;
         searchTree.getNearestNeighbour(gLocations[i].lat(), gLocations[i].lon(), Y, X);
         gYi[i] = Y;
         gXi[i] = X;
         gLafs[i] = lafs[Y][X];

         // Check if the elevation of the station roughly matches the reference grid elevation
         bool hasValidElev = Util::isValid(gLocations[i].elev());
         if(Util::isValid(mMaxElevDiff) && hasValidElev) {
            float elevDiff = abs(gLocations[i].elev() - elevs[Y][X]);
            hasValidElev = elevDiff < mMaxElevDiff;
         }
         if(hasValidElev) {
            for(int y = std::max(0, Y - radius); y < std::min(nY, Y + radius); y++) {
               for(int x = std::max(0, X - radius); x < std::min(nX, X + radius); x++) {
                  if(mSort || gLocIndices[y][x].size() < mMaxLocations) {
                     gLocIndices[y][x].push_back(i);
                     count ++;
                  }
               }
            }
         }
         else {
            std::stringstream ss;
            ss << "Removing station because elevation (" << gLocations[i].elev() << " m) is too far from grid (" << elevs[Y][X] << " m)";
            Util::warning(ss.str());
         }
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
      if(mNumVariable != "") {
         num = iFile.getField(mNumVariable, t);
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               for(int e = 0; e < nEns; e++) {
                  (*num)(y, x, e) = 0;
               }
            }
         }
      }

      // Compute Y
      vec2 gY(gS);
      std::vector<float> gYhat(gS);
      for(int i = 0; i < gS; i++) {
         gY[i].resize(nEns, 0);
         float elevCorr = 0;
         if(Util::isValid(mElevGradient)) {
            float nnElev = elevs[gYi[i]][gXi[i]];
            assert(Util::isValid(nnElev));
            assert(Util::isValid(gLocations[i].elev()));
            float elevDiff = gLocations[i].elev() - nnElev;
            elevCorr = mElevGradient * elevDiff;
         }
         float total = 0;
         int count = 0;
         for(int e = 0; e < nEns; e++) {
            float value = (*field)(gYi[i], gXi[i], e);
            if(Util::isValid(value)) {
               value += elevCorr;
               gY[i][e] = value;
               total += value;
               count++;
            }
         }
         float mean = total / count;
         for(int e = 0; e < nEns; e++) {
            float value = gY[i][e];
            if(Util::isValid(value)) {
               gY[i][e] -= mean;
            }
         }
         gYhat[i] = mean;
         if(useBias) {
            float currBias = (*bias)(gYi[i], gXi[i], 0);
            gYhat[i] -= currBias;
         }
      }

      #pragma omp parallel for
      for(int x = 0; x < nX; x++) {
         for(int y = 0; y < nY; y++) {
            float lat = lats[y][x];
            float lon = lons[y][x];
            float elev = elevs[y][x];
            float laf = lafs[y][x];

            //
            // Create list of locations for this gridpoint
            //
            std::vector<int> lLocIndices0 = gLocIndices[y][x];
            std::vector<int> lLocIndices;
            lLocIndices.reserve(lLocIndices0.size());
            std::vector<std::pair<float,int> > lRhos0;
            lRhos0.reserve(lLocIndices0.size());
            for(int i = 0; i < lLocIndices0.size(); i++) {
               int index = lLocIndices0[i];
               float hdist = Util::getDistance(gLocations[index].lat(), gLocations[index].lon(), lat, lon, true);
               float vdist = Util::MV;
               if(Util::isValid(gLocations[index].elev() && Util::isValid(elev)))
                  vdist = gLocations[index].elev() - elev;
               float lafdist = 0;
               if(Util::isValid(gLafs[index]) && Util::isValid(laf))
                  lafdist = gLafs[index] - laf;
               float rho = calcRho(hdist, vdist, lafdist);
               int X = gXi[index];
               int Y = gYi[index];
               // Only include observations that are within the domain
               if(X > 0 && X < lats[0].size()-1 && Y > 0 && Y < lats.size()-1) {
                  if(rho > mMinRho) {
                     lRhos0.push_back(std::pair<float,int>(rho, i));
                  }
               }
            }
            arma::vec lRhos;
            if(lRhos0.size() > mMaxLocations) {
               // If sorting is enabled and we have too many locations, then only keep the best ones based on rho.
               // Otherwise, just use the last locations added
               lRhos = arma::vec(mMaxLocations);
               if(mSort) {
                  std::sort(lRhos0.begin(), lRhos0.end(), Util::sort_pair_first<float,int>());
               }
               for(int i = 0; i < mMaxLocations; i++) {
                  // The best values start at the end of the array
                  int index = lRhos0[lRhos0.size() - 1 - i].second;
                  lLocIndices.push_back(lLocIndices0[index]);
                  lRhos(i) = lRhos0[lRhos0.size() - 1 - i].first;
               }
            }
            else {
               lRhos = arma::vec(lRhos0.size());
               for(int i = 0; i < lRhos0.size(); i++) {
                  int index = lRhos0[i].second;
                  lLocIndices.push_back(lLocIndices0[index]);
                  lRhos(i) = lRhos0[i].first;
               }
            }

            int lS = lLocIndices.size();

            if(lS == 0 || lS < mMinObs) {
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
               std::vector<float> lObs(lS, Util::MV);
               for(int i = 0; i < lS; i++) {
                  int index = lLocIndices[i];
                  lObs[i] = gObs[index];
               }
               float value = Util::calculateStat(lObs, Util::StatTypeQuantile, 0.5);
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

            vectype lObs(lS);
            vectype lElevs(lS);
            vectype lLafs(lS);
            for(int i = 0; i < lLocIndices.size(); i++) {
               int index = lLocIndices[i];
               lObs[i] = gObs[index];
               lElevs[i] = gElevs[index];
               lLafs[i] = gLafs[index];
            }

            // Compute Y (model at obs-locations)
            mattype lY(lS, nValidEns);
            vectype lYhat(lS);

            for(int i = 0; i < lS; i++) {
               // Use the nearest neighbour for this location
               int index = lLocIndices[i];
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  lY(i, e) = gY[index][ei];
               }
               lYhat(i) = gYhat[index];

            }

            // 
            // Revert to static structure function when there is not enough ensemble information
            // 
            if(!mUseEns || nValidEns < mMinValidEns) {
               // Current grid-point to station error covariance matrix
               mattype lG(1, lS, arma::fill::zeros);
               // Station to station error covariance matrix
               mattype lP(lS, lS, arma::fill::zeros);
               // Station variance
               mattype lR(lS, lS, arma::fill::zeros);
               for(int i = 0; i < lS; i++) {
                  int index = lLocIndices[i];
                  lR(i, i) = gCi[index];
                  float hdist = Util::getDistance(gLocations[index].lat(), gLocations[index].lon(), lat, lon, true);
                  float vdist = Util::MV;
                  if(Util::isValid(gLocations[index].elev() && Util::isValid(elev)))
                     vdist = gLocations[index].elev() - elev;
                  float lafdist = 0;
                  if(Util::isValid(gLafs[index]) && Util::isValid(laf))
                     lafdist = gLafs[index] - laf;
                  float rho = calcRho(hdist, vdist, lafdist);
                  lG(0, i) = rho;
                  for(int j = 0; j < lS; j++) {
                     int index_j = lLocIndices[j];
                     float hdist = Util::getDistance(gLocations[index].lat(), gLocations[index].lon(), gLocations[index_j].lat(), gLocations[index_j].lon(), true);
                     float vdist = Util::MV;
                     if(Util::isValid(gLocations[index].elev() && Util::isValid(gLocations[index_j].elev())))
                        vdist = gLocations[index].elev() - gLocations[index_j].elev();
                           float lafdist = 0;
                     if(Util::isValid(gLafs[index]) && Util::isValid(laf))
                        lafdist = gLafs[index] - gLafs[index_j];

                     lP(i, j) = calcRho(hdist, vdist, lafdist);
                  }
               }
               mattype lGSR;
               if(useBias)
                  lGSR = lG * arma::inv(lP + 1 / (1 + mGamma) * mEpsilon * mEpsilon * lR);
               else
                  lGSR = lG * arma::inv(lP + mEpsilon * mEpsilon * lR);

               // This should loop over nValidEns. And use ei.
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  if (Util::isValid((*field)(y, x, ei))) {
                     vectype currFcst = lY.col(e) + lYhat;
                     vectype dx = lGSR * (lObs - currFcst);
                     if(x == mX && y == mY) {
                        std::cout << "Lat: " << lat << std::endl;
                        std::cout << "Lon: " << lon << " " << lat << " " << std::endl;
                        std::cout << "Elev: " << elev << std::endl;
                        std::cout << "Elev: " << elev << std::endl;
                        std::cout << "P:" << std::endl;
                        print_matrix<mattype>(lP);
                        std::cout << "R:" << std::endl;
                        print_matrix<mattype>(lR);
                        std::cout << "GSR:" << std::endl;
                        print_matrix<mattype>(lGSR);
                        std::cout << "Obs:" << std::endl;
                        print_matrix<mattype>(lObs);
                        std::cout << "Current forecast: " << std::endl;
                        print_matrix<mattype>(currFcst);
                        std::cout << "Increment" << std::endl;
                        print_matrix<mattype>(lObs - currFcst);
                        std::cout << "Yhat" << std::endl;
                        print_matrix<mattype>(lYhat);
                        std::cout << "dx: " << dx[0] << std::endl;
                     }
                     (*output)(y, x, ei) = (*field)(y, x, ei) + dx[0];
                  }
               }

               // Update bias
               if(useBias) {
                  float biasTotal = 0;
                  (*newbias)(y, x, 0) = (*bias)(y, x, 0) - mGamma / (1 + mGamma) * biasTotal;
               }
            }
            else {
               // Compute Rinv
               mattype Rinv(lS, lS, arma::fill::zeros);
               for(int i = 0; i < lS; i++) {
                  int index = lLocIndices[i];
                  float r = sigma * sigma * gCi[index];
                  if(mUseRho) {
                     float rho = lRhos[i];
                     Rinv(i, i) = 1 / r * rho;
                     if(x == mX && y == mY) {
                        std::cout << "R(" << i << ") " << Rinv(i, i) << std::endl;
                     }
                  }
                  else {
                     Rinv(i, i) = 1 / r;
                  }
               }

               // Compute C matrix
               // k x gS * gS x gS
               mattype C(nValidEns, lS);
               C = lY.t() * Rinv;

               mattype Pinv(nValidEns, nValidEns);
               float currDelta = 1;
               if(Util::isValid(mDelta))
                  currDelta = mDelta;
               else if(useDelta) {
                  currDelta = (*delta)(y, x, 0);
               }
               float diag = 1 / currDelta * (nValidEns - 1);
               if(useBias)
                  diag = 1 / currDelta / (1 + mGamma) * (nValidEns - 1);

               Pinv = C * lY + diag * arma::eye<mattype>(nValidEns, nValidEns);
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
                  print_matrix<mattype>(lY);
                  std::cout << "lObs:" << std::endl;
                  print_matrix<mattype>(lObs);
                  std::cout << "Yhat" << std::endl;
                  print_matrix<mattype>(lYhat);
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
               mattype PC(nValidEns, lS);
               PC = P * C;

               // Compute w
               vectype w(nValidEns);
               if(mDiagnose)
                  w = PC * (arma::ones<vectype>(lS));
               else
                  w = PC * (lObs - lYhat);

               // Add w to W
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
               float ensMean = total / count;
               for(int e = 0; e < nValidEns; e++) {
                  X(e) -= ensMean;
               }

               // Write debugging information
               if(x == mX && y == mY) {
                  std::cout << "Lat: " << lat << std::endl;
                  std::cout << "Lon: " << lon << " " << lat << " " << std::endl;
                  std::cout << "Elev: " << elev << std::endl;
                  std::cout << "Laf: " << laf << std::endl;
                  std::cout << "Num obs: " << lS << std::endl;
                  std::cout << "Num ens: " << nValidEns << std::endl;
                  std::cout << "rhos" << std::endl;
                  print_matrix<mattype>(lRhos);
                  std::cout << "P" << std::endl;
                  print_matrix<mattype>(P);
                  std::cout << "C" << std::endl;
                  print_matrix<mattype>(C);
                  std::cout << "C * lY" << std::endl;
                  print_matrix<mattype>(C * lY);
                  std::cout << "PC" << std::endl;
                  print_matrix<mattype>(PC);
                  std::cout << "W" << std::endl;
                  print_matrix<mattype>(W);
                  std::cout << "w" << std::endl;
                  print_matrix<mattype>(w);
                  std::cout << "Y:" << std::endl;
                  print_matrix<mattype>(lY);
                  std::cout << "Yhat" << std::endl;
                  print_matrix<mattype>(lYhat);
                  std::cout << "lObs" << std::endl;
                  print_matrix<mattype>(lObs);
                  std::cout << "lObs - Yhat" << std::endl;
                  print_matrix<mattype>(lObs - lYhat);
                  std::cout << "X" << std::endl;
                  print_matrix<mattype>(X);
                  std::cout << "elevs" << std::endl;
                  print_matrix<mattype>(lElevs);
                  std::cout << "lafs" << std::endl;
                  print_matrix<mattype>(lLafs);
                  std::cout << "Analysis increment:" << std::endl;
                  print_matrix<mattype>(X.t() * W);
                  std::cout << "My: " << arma::mean(arma::dot(lObs - lYhat, lRhos) / lS) << std::endl;
               }

               // Compute analysis
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  float total = 0;
                  for(int k = 0; k < nValidEns; k++) {
                     total += X(k) * W(k, e);
                  }

                  float currIncrement = total;
                  if(mUseMeanBias)
                     currIncrement = arma::mean(lObs - lYhat);

                  if(mSaveDiff)
                     (*output)(y, x, ei) = currIncrement;
                  else {
                     float raw = ensMean;
                     if(useBias) {
                        raw -= (*bias)(y, x, 0);
                     }
                     if(mUseBug) {
                        raw = (*field)(y, x, ei);
                     }
                     if(!mExtrapolate) {
                        // Don't allow a final increment that is larger than any increment
                        // at station points
                        float maxInc = arma::max(lObs - (lY[e] + lYhat));
                        float minInc = arma::min(lObs - (lY[e] + lYhat));
                        if(x == mX && y == mY) {
                           std::cout << "Increments: " << maxInc << " " << minInc << " " << currIncrement << std::endl;
                        }

                        // The increment for this member. currIncrement is the increment relative to
                        // ensemble mean
                        float memberIncrement = currIncrement - X(e);
                        // Adjust increment if it gives a member increment that is outside the range
                        // of the observation increments
                        if(x == mX && y == mY) {
                           std::cout << "Analysis increment: " << memberIncrement << " " << ensMean << " " << currIncrement << " " << X(e) << std::endl;
                        }
                        if(maxInc > 0 && memberIncrement > maxInc) {
                           currIncrement = maxInc + X(e);
                        }
                        else if(maxInc < 0 && memberIncrement > 0) {
                           currIncrement = 0 + X(e);
                        }
                        else if(minInc < 0 && memberIncrement < minInc) {
                           currIncrement = minInc + X(e);
                        }
                        else if(minInc > 0 && memberIncrement < 0) {
                           currIncrement = 0 + X(e);
                        }
                        if(x == mX && y == mY) {
                           std::cout << "Final increment: " << currIncrement << " " << currIncrement - X(e) << std::endl;
                        }
                     }
                     (*output)(y, x, ei) = ensMean + currIncrement;
                  }

                  if(mNumVariable != "") {
                     (*num)(y, x, ei) = lS;
                  }
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
               float trace = arma::trace(lY * lY.t());
               float numerator = mSigma * mSigma / mEpsilon / mEpsilon;
               float denomenator = 1.0 / lS / (nValidEns - 1) * trace;
               float currDeltaEvidence = numerator / denomenator;
               float weightOld = deltaVar;
               float weightNew = mNewDeltaVar;
               (*newdelta)(y, x, 0) = ((*delta)(y, x, 0) * weightNew + currDeltaEvidence * weightOld) / (weightOld + weightNew);
               */
            }
         }
      }
      iFile.addField(output, mVariable, t);
      if(mNumVariable != "")
         iFile.addField(num, Variable(mNumVariable), t);
      if(useBias) {
         iFile.addField(newbias, Variable(mBiasVariable), t);
      }
      if(useDelta) {
         float oldDelta = (*delta)(0, 0, 0);
         float value = calcDelta(oldDelta, gY);
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               (*newdelta)(y, x, 0) = value;
            }
         }

         iFile.addField(newdelta, Variable(mDeltaVariable), t);
      }
   }
   return true;
}

float CalibratorOi::calcDelta(float iOldDelta, const vec2& iY) const {
   float deltaVar = mC - 1;
   float trace = 0;
   float numValidS = 0;
   int S = iY.size();
   int nEns = iY[0].size();
   for(int s = 0; s < S; s++) {
      // Compute value in coord s,s
      float value = 0;
      int count = 0;
      for(int e = 0; e < nEns; e++) {
         if(Util::isValid(iY[s][e])) {
            value += iY[s][e] * iY[s][e];
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
   return (iOldDelta * weightNew + currDeltaEvidence * weightOld) / (weightOld + weightNew);
}

float CalibratorOi::calcRho(float iHDist, float iVDist, float iLDist) const {
   float h = (iHDist/mHLength);
   float rho = exp(-0.5 * h * h);
   if(Util::isValid(mVLength)) {
      if(!Util::isValid(iVDist)) {
         rho = 0;
      }
      else {
         float v = (iVDist/mVLength);
         rho *= exp(-0.5 * v * v);
      }
   }
   if(Util::isValid(mWMin)) {
      rho *= 1 - (1 - mWMin) * std::abs(iLDist);
   }
   return rho;
}

std::string CalibratorOi::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c oi","Spreads bias in space by using kriging. A parameter file is required, which must have one column with the bias.")<< std::endl;
   return ss.str();
}
