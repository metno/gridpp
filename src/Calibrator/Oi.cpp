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
      mMu(0.25),
      mSort(true),
      mBiasVariable(""),
      mMaxLocations(20),
      mGamma(0.9) {
   iOptions.getValue("bias", mBiasVariable);
   iOptions.getValue("d", mHLength);
   iOptions.getValue("h", mVLength);
   iOptions.getValue("maxLocations", mMaxLocations);
   iOptions.getValue("sort", mSort);
}

template<class Matrix>
void print_matrix(Matrix matrix) {
       matrix.print(std::cout);
}
template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);

bool CalibratorOi::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

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

   // Create nearest neighbour tree
   std::vector<float> ci(obsLocations.size());
   std::vector<float> obs(obsLocations.size());
   for(int i = 0; i < obsLocations.size(); i++) {
      Parameters parameters = iParameterFile->getParameters(0, obsLocations[i]);
      obs[i] = parameters[0];
      ci[i] = parameters[1];
   }
   std::cout << "Loaded parameters" << std::endl;
   float sigma = 2;
   int S = obsLocations.size();
   std::vector<std::vector<std::vector<int> > > obsIndices; // Y, X, obs indices
   std::vector<float> obsY(S);
   std::vector<float> obsX(S);
   obsIndices.resize(nY);
   for(int y = 0; y < nY; y++) {
      obsIndices[y].resize(nX);
      for(int x = 0; x < nX; x++) {
   /*
         obsIndices[y][x].reserve(mMaxLocations);
  */
      }
   }

   int radius = 3.64 * mHLength / 1000; // 2*mHLength / 1000;

   double time_s = Util::clock();
   KDTree searchTree(iFile.getLats(), iFile.getLons());
   int count = 0;
   std::cout << "Setting up locations radius=" << radius << std::endl;
   for(int i = 0; i < S; i++) {
      if(i % 1000 == 0)
         std::cout << i << std::endl;
      // if(Util::isValid(obs[i]) && obsLocations[i].lat() < 62 && obsLocations[i].lat() > 58 && obsLocations[i].lon() > 8 && obsLocations[i].lon() < 12) {
      if(Util::isValid(obs[i])) {
         int Y, X;
         searchTree.getNearestNeighbour(obsLocations[i].lat(), obsLocations[i].lon(), Y, X);
         for(int y = std::max(0, Y - radius); y < std::min(nY, Y + radius); y++) {
            for(int x = std::max(0, X - radius); x < std::min(nX, X + radius); x++) {
               // if(obsIndices[y][x].size() == 0)
               //    obsIndices[y][x].resize(mMaxLocations);
               if(1 || obsIndices[y][x].size() < mMaxLocations) {
                  obsIndices[y][x].push_back(i);
                  count++;
               }
            }
         }
         obsY[i] = Y;
         obsX[i] = X;
      }
   }
   double time_e = Util::clock();
   std::cout << "Finding locations " << time_e - time_s << " " << count << std::endl;

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      FieldPtr output = iFile.getEmptyField();
      FieldPtr bias;
      if(iFile.hasVariable(mBiasVariable))
         bias = iFile.getField(mBiasVariable, t);
      FieldPtr accum = iFile.getEmptyField(0);
      FieldPtr weights = iFile.getEmptyField(0);

#if 0
      #pragma omp parallel for
      for(int yx = 0; yx < nY*nX; yx++) {
         int y = yx / nX;
         int x = yx % nX;
         {
#else
      for(int y = 0; y < nY; y++) {
         double time_s = Util::clock();
         float totalLocations = 0;
         for(int x = 0; x < nX; x++) {
#endif
            float lat = lats[y][x];
            float lon = lons[y][x];
            float elev = elevs[y][x];
            std::stringstream ss;
            std::vector<int> useLocations0 = obsIndices[y][x];
            std::vector<int> useLocations;
            useLocations.reserve(useLocations0.size());
            std::vector<std::pair<float,int> > rhos(useLocations0.size());
            for(int i = 0; i < useLocations0.size(); i++) {
               int index = useLocations0[i];
               float hdist = Util::getDistance(obsLocations[index].lat(), obsLocations[index].lon(), lat, lon, true);
               float vdist = obsLocations[index].elev() - elev;
               float rho = calcRho(hdist, vdist);
               int X = obsX[index];
               int Y = obsY[index];
               if (x == 650 && y == 770) {
                  std::cout << "### " << i << " " << vdist << " " << rho << std::endl;
               }
               if(X > 0 && X < lats[0].size() && Y > 0 && Y < lats.size() and useLocations.size() <= mMaxLocations) {
                  if(rho > 0.0013) {
                     rhos[i] = std::pair<float,int>(rho, i);
                     // useLocations.push_back(useLocations0[i]);
                  }
               }
            }
            if(useLocations0.size() > mMaxLocations) {
               if(mSort)
                  std::sort(rhos.begin(), rhos.end(), Util::sort_pair_first<float,int>());
               for(int i = 0; i < mMaxLocations; i++) {
                  int index = rhos[rhos.size() - 1 - i].second;
                  useLocations.push_back(useLocations0[index]);
               }
            }
            else {
               for(int i = 0; i < rhos.size(); i++) {
                  int index = rhos[i].second;
                  useLocations.push_back(useLocations0[index]);
               }
            }

            if (x == 650 && y == 770) {
               std::cout << "### " << useLocations0.size() << " " << useLocations.size() << std::endl;

            }
            assert(useLocations.size() < 10000);
            int nObs = useLocations.size();
            // if(nObs <= 3 || x < 10 || y < 10)
            bool allSame = true;
            for(int  i = 1; i < nObs; i++) {
               allSame = allSame && obsX[useLocations[i]] == obsX[useLocations[i-1]] && obsY[useLocations[i]] == obsY[useLocations[i-1]];
            }
            if(nObs <= 3 || allSame) {
               for(int e = 0; e < nEns; e++) {
                  (*output)(y, x, e) = (*field)(y, x, e); // Util::MV;
               }
               continue;
            }
            // std::cout << useLocations.size() << std::endl;

            // Compute Rinv
            arma::vec rinv(nObs);
            for(int i = 0; i < useLocations.size(); i++) {
               int index = useLocations[i];
               float dist = Util::getDistance(obsLocations[index].lat(), obsLocations[index].lon(), lat, lon, true);
               float vdist = obsLocations[index].elev() - elev;
               float rho = calcRho(dist, vdist);
               float r = sigma * ci[index];
               // rvec(i) = r;
               rinv(i) = 1 / r * rho;
            }
            arma::mat Rinv = arma::diagmat(rinv);

            // Compute Y (model at obs-locations)
            arma::mat Y(nObs, nEns);
            arma::vec Yhat(nObs);
            for(int i = 0; i < nObs; i++) {
               // Use the nearest neighbour for this location
               float total = 0;
               int count = 0;
               for(int e = 0; e < nEns; e++) {
                  float value = (*field)(obsY[i], obsX[i], e);
                  if(Util::isValid(value)) {
                     Y(i, e) = value;
                     total += value;
                     count++;
                  }
               }
               // assert(count > 0);
               float mean = total / count;
               Yhat(i) = mean;
               for(int e = 0; e < nEns; e++) {
                  Y(i, e) -= mean;
               }
            }

            // Compute C matrix
            // k x S * S x S
            arma::mat C(nEns, nObs);
            C = Y.t() * Rinv;
            /*
            for(int e = 0; e < nEns; e++) {
               for(int i = 0; i < nObs; i++) {
                  C(e, i) += Y(i, e) * Rinv(i, i);
               }
            }
            */

            arma::mat Pinv(nEns, nEns);
            float delta = 1;
            float diag = 1/delta / (1 + mGamma) * (nEns - 1);
            Pinv = C * Y + diag * arma::eye(nEns, nEns);
            /*
            for(int e = 0; e < nEns; e++) {
               for(int e2 = 0; e2 < nEns; e2++) {
                  for(int k = 0; k < useLocations.size(); k++) {
                     Pinv(e, e2) += C(e, k) * Y(k, e2);
                  }
                  if(e == e2)
                     Pinv(e,e2) += diag;
               }
            }
            */
            // arma::Mat<float> P = arma::inv(Pinv);
            float cond = arma::rcond(Pinv);
            // std::cout << "Condition number: " << cond << std::endl;
            if(cond <= 0) {
               for(int e = 0; e < nEns; e++) {
                  (*output)(y, x, e) = (*field)(y, x, e); // Util::MV;
               }
               // print_matrix<arma::mat>(Pinv);
               continue;
            }
            arma::mat P = arma::inv(Pinv);
            bool issingular = false;
            for(int e = 0; e < nEns; e++) {
               if(P(e, e) == 0) {
                  issingular = true;
                  break;
               }
            }
            if(issingular) {
               std::cout << "Singular" << std::endl;
               // for(int e = 0; e < nEns; e++) {
               //    (*output)(y, x, e) = Util::MV;
               // }
               // continue;
            }
            arma::mat W = arma::real(arma::sqrtmat((nEns - 1) * P));
            if(W.n_rows == 0) {
               for(int e = 0; e < nEns; e++) {
                  (*output)(y, x, e) = (*field)(y, x, e); // Util::MV;
               }
               continue;
            }

            /*
            // Compute P matrix
            vec2 P(nEns);
            for(int e = 0; e < nEns; e++) {
               P[e].resize(nEns, 0);
               float diag = 1/delta / (1 + mGamma) * (nEns - 1);
               for(int e2 = 0; e2 < nEns; e2++) {
                  // TODO account for missing ensemble members
                  for(int k = 0; k < useLocations.size(); k++) {
                     P[e][e2] += C[e][k] * Y[k][e2];
                  }
                  if(e == e2)
                     P[e][e2] += diag;
               }
            }
            P = Util::inverse(P);
            */

            // Compute PC
            arma::mat PC(nEns, nObs);
            for(int e = 0; e < nEns; e++) {
               for(int i = 0; i < nObs; i++) {
                  for(int k = 0; k < nEns; k++) {
                     // PC[e][i] += P[e][k] * C[k][i];
                     PC(e, i) += P(e, k) * C(k, i);
                  }
               }
            }

            // Compute w
            arma::vec w(nEns);
            for(int e = 0; e < nEns; e++) {
               for(int i = 0; i < useLocations.size(); i++) {
                  int index = useLocations[i];
                  float currObs = obs[index];
                  float dy = currObs - Yhat(i);
                  if(dy > 3)
                     dy = 3;
                  if(dy < -3)
                     dy = -3;
                  assert(abs(dy) < 10);
                  w(e) += PC(e, i) * dy;
                  //w[e] += currObs / useLocations.size();
               }
            }


            // Compute X (perturbations about model mean)
            arma::vec X(nEns);
            float total = 0;
            int count = 0;
            for(int e = 0; e < nEns; e++) {
               float value = (*field)(y, x, e);
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
            // assert(count > 0);
            float total2 = 0;
            for(int e = 0; e < nEns; e++) {
               X(e) -= mean;
               total2 += X(e);
            }
            // if(total2 != 0) {
            //    print_matrix<arma::mat>(X);
            //    std::cout << total2 << std::endl;
            // }

            // assert(total2 == 0);
            // assert(X(0) != 0);

            if(0 && nObs == 4) {
               print_matrix<arma::mat>(P);
               print_matrix<arma::mat>(W);
               print_matrix<arma::mat>(X);
               print_matrix<arma::mat>(W);
               std::cout << "Analaysis:" << std::endl;
               print_matrix<arma::mat>(X.t() * W);
            }
            // Compute analysis
            for(int e = 0; e < nEns; e++) {
               float total = 0;
               for(int k = 0; k < nEns; k++) {
                  total += X(k) * W(k, e);
                  if(0 && nObs == 4)
                     std::cout << "# " << total << std::endl;

               }
               // assert(total < 100);
               // assert(total > -100);
               (*output)(y, x, e) = (*field)(y, x, e) + total;
               if(total == 0) {
                  std::cout << "P" << std::endl;
                  print_matrix<arma::mat>(P);
                  std::cout << "W" << std::endl;
                  print_matrix<arma::mat>(W);
                  std::cout << "Y:" << std::endl;
                  print_matrix<arma::mat>(Y);
                  std::cout << "Yhat" << std::endl;
                  print_matrix<arma::mat>(Yhat);
                  std::cout << "X" << std::endl;
                  print_matrix<arma::mat>(X);
                  std::cout << "Analaysis:" << std::endl;
                  print_matrix<arma::mat>(X.t() * W);
                  abort();
               }
               if(0 && nObs == 4)
                  std::cout << " " << total << std::endl;
               // (*output)(y, x, e) = useLocations.size();
               // for(int k = 0; k < nEns; k++) {
               //    (*field)(y, x, e) += P[e][k];
               // }
               // (*field)(y, x, e) = w[e];
            }

            /*
            const Location currLocation(lats[y][x], lons[y][x]);
            std::vector<int> useLocations;
            useLocations.reserve(obsLocations.size());

            // Compute observation error covariance matrix

            for(int i = 0; i < obsLocations.size(); i++) {
               float dist = Util::getDistance(obsLocations[i].lat(), obsLocations[i].lon(), lats[y][x], lons[y][x], true);
               if(useLocations.size() > 20)
                  break;
               if(dist < mHLength) {
                  useLocations.push_back(i);
               }
            }
            */
            totalLocations += useLocations.size()*useLocations.size();

            /*
            int S = useLocations.size();
            ss << "Computing " << y << " " << x << " " << S;
            Util::info(ss.str());
            vec2 Rmatrix(S);
            // vec2 Rmatrix_inv(S);
            for(int i = 0; i < S; i++) {
               Rmatrix[i].resize(S, 0);
               Rmatrix[i][i] = sigma * ci[i];
               //Rmatrix_inv[i].resize(S, 0);
               //Rmatrix_inv[i][i] = 1 / (sigma * ci[i]);
            }
            vec2 Rmatrix_inv = Util::inverse(Rmatrix);

            // Location nearestLocation(Util::MV, Util::MV);
            // iParameterFile->getNearestLocation(t, currLocation, nearestLocation);
            const Parameters parameters = iParameterFile->getParameters(t, currLocation);
            for(int e = 0; e < nEns; e++) {
               (*field)(y, x, e) = parameters[0];
            }
  */
         }
         double time_e = Util::clock();
         // std::cout << y << " " << time_e - time_s << "s expected total time: " << (time_e - time_s) * nX << std::endl; // << " " << totalLocations / nX << " " << (time_e - time_s) / (totalLocations / nX)<< std::endl;
      }
      iFile.addField(output, mVariable, t);
   }
   return true;
}

float CalibratorOi::calcRho(float iHDist, float iVDist) const {
   // float h = (mHLength*mHLength - iHDist*iHDist) / (mHLength*mHLength + iHDist*iHDist);
   // float v = (mVLength*mVLength - iVDist*iVDist) / (mVLength*mVLength + iVDist*iVDist);
   float h = (iHDist/mHLength);
   float v = (iVDist/mVLength);
   float rho = exp(-0.5 * h * h) * exp(-0.5 * v * v);
   return rho;
}
float CalibratorOi::calcCovar(const Location& loc1, const Location& loc2) const {
   float weight = Util::MV;
   bool isValidLoc1 = Util::isValid(loc1.lat()) && Util::isValid(loc1.lon()) && Util::isValid(loc1.elev());
   bool isValidLoc2 = Util::isValid(loc2.lat()) && Util::isValid(loc2.lon()) && Util::isValid(loc2.elev());
   float useApproxDistance = true;
   if(isValidLoc1 && isValidLoc2) {
      float horizDist = Util::getDistance(loc1.lat(), loc1.lon(), loc2.lat(), loc2.lon(), useApproxDistance);
      float vertDist  = fabs(loc1.elev() - loc2.elev());
      if(horizDist >= mHLength || (Util::isValid(mVLength) && vertDist >= mVLength)) {
         // Outside radius of influence
         return 0;
      }
      float horizWeight = (mHLength*mHLength - horizDist*horizDist) /
                          (mHLength*mHLength + horizDist*horizDist);
      float vertWeight = 1;
      if(Util::isValid(mVLength)) {
         vertWeight  = (mVLength*mVLength - vertDist*vertDist) /
                             (mVLength*mVLength + vertDist*vertDist);
      }
      weight = horizWeight * vertWeight;
   }
   return weight;
}

std::string CalibratorOi::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c oi","Spreads bias in space by using kriging. A parameter file is required, which must have one column with the bias.")<< std::endl;
   return ss.str();
}
