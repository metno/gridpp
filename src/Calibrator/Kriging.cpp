#include "Kriging.h"
#include "../Util.h"
#include "../Parameters.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <math.h>

CalibratorKriging::CalibratorKriging(Variable::Type iVariable, const Options& iOptions):
      Calibrator(iOptions),
      mVariable(iVariable),
      mEfoldDist(30000),
      mMaxElevDiff(100),
      mAuxVariable(Variable::None),
      mLowerThreshold(Util::MV),
      mUpperThreshold(Util::MV),
      mRadius(30000),
      mOperator(Util::OperatorAdd),
      // Use an approximation when calculating distances between points?
      mUseApproxDistance(true),
      mKrigingType(TypeCressman) {
   iOptions.getValue("efoldDist", mEfoldDist);
   iOptions.getValue("radius", mRadius);
   iOptions.getValue("maxElevDiff", mMaxElevDiff);
   std::string parfilename;
   if(mEfoldDist< 0) {
      Util::error("CalibratorKriging: 'efoldDist' must be >= 0");
   }
   if(mMaxElevDiff < 0) {
      Util::error("CalibratorKriging: 'maxElevDiff' must be >= 0");
   }
   if(mRadius < 0) {
      Util::error("CalibratorKriging: 'radius' must be >= 0");
   }

   std::string type;
   if(iOptions.getValue("type", type)) {
      if(type == "cressman") {
         mKrigingType = TypeCressman;
      }
      else if(type == "barnes") {
         mKrigingType = TypeBarnes;
      }
      else {
         Util::error("CalibratorKriging: 'type' not recognized");
      }
   }

   std::string op;
   if(iOptions.getValue("operator", op)) {
      if(op == "add") {
         mOperator = Util::OperatorAdd;
      }
      else if(op == "subtract") {
         mOperator = Util::OperatorSubtract;
      }
      else if(op == "multiply") {
         mOperator = Util::OperatorMultiply;
      }
      else if(op == "divide") {
         mOperator = Util::OperatorDivide;
      }
      else {
         Util::error("CalibratorKriging: 'operator' not recognized");
      }
   }

   std::string auxVariable;
   if(iOptions.getValue("auxVariable", auxVariable)) {
      mAuxVariable = Variable::getType(auxVariable);
      std::vector<float> range;
      if(iOptions.getValues("range", range)) {
         if(range.size() != 2) {
            Util::error("CalibratorKriging: 'range' must be of the form lower,upper");
         }
         mLowerThreshold = range[0];
         mUpperThreshold = range[1];
      }
      else {
         Util::error("CalibratorKriging: 'range' required if using 'auxVariable'.");
      }
      if(mLowerThreshold > mUpperThreshold) {
         Util::error("CalibratorKriging: the lower value must be less than upper value in 'range'");
      }
   }
}

bool CalibratorKriging::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   if(!iParameterFile->isLocationDependent()) {
      std::stringstream ss;
      ss << "Kriging requires a parameter file with spatial information";
      Util::error(ss.str());
   }

   // Precompute if a gridpoint has an observation location (for any point in time) within the
   // radius of influence. This saves time for gridpoints without observations close by, since this
   // does not need to be determined for each forecast time.
   std::vector<Location> pointLocations = iParameterFile->getLocations();
   std::vector<std::vector<bool> > hasObs;
   hasObs.resize(nLat);
   int total = 0;
   #pragma omp parallel for  reduction(+:total)
   for(int i = 0; i < nLat; i++) {
      hasObs[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         hasObs[i][j] = false;
         float lat = lats[i][j];
         float lon = lons[i][j];
         float elev = elevs[i][j];
         for(int k = 0; k < pointLocations.size(); k++) {
            Location loc = pointLocations[k];
            if(Util::getDistance(loc.lat(), loc.lon(), lat, lon, mUseApproxDistance) <= mRadius && fabs(loc.elev() - elev) <= mMaxElevDiff) {
               hasObs[i][j] = true;
               total++;
               break;
            }
         }
      }
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      FieldPtr accum = iFile.getEmptyField(0);
      FieldPtr weights = iFile.getEmptyField(0);
      FieldPtr auxField;
      if(mAuxVariable != Variable::None)
         auxField = iFile.getField(mAuxVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            if(!hasObs[i][j]) {
               continue;
            }
            float lat = lats[i][j];
            float lon = lons[i][j];
            float elev = elevs[i][j];

            // Find stations within radius
            std::vector<Location> validLocations;
            std::vector<float> bias;
            for(int k = 0; k < pointLocations.size(); k++) {
               Location loc = pointLocations[k];
               if(Util::getDistance(loc.lat(), loc.lon(), lat, lon, mUseApproxDistance) <= mRadius && fabs(loc.elev() - elev) <= mMaxElevDiff) {
                  Parameters parameters = iParameterFile->getParameters(t, loc);
                  if(parameters.size() > 0) {
                     validLocations.push_back(loc);
                     float currBias = parameters[0];
                     if(mOperator == Util::OperatorMultiply) {
                        currBias = currBias - 1;
                     }
                     else if(mOperator == Util::OperatorDivide) {
                        currBias = currBias - 1;
                     }
                     bias.push_back(currBias);
                  }
               }
            }
            int N = validLocations.size();
            if(N > 0) {
               // Set up matrix system: 
               // matrix * weights = S
               // matrix: The obs-to-obs weight (NxN)
               // S: The obs-to-current_grid_point weight (Nx1)
               // bias: The bias at each obs location (Nx1)
               // weights = (matrix)^-1 * S
               // gridpoint_bias = weights' * bias (scalar)
               vec2 matrix;
               matrix.resize(N);
               for(int ii = 0; ii < N; ii++) {
                  matrix[ii].resize(N,0);
               }
               std::vector<float> S;
               S.resize(N);
               for(int ii = 0; ii < N; ii++) {
                  Location iloc = validLocations[ii];

                  S[ii] = calcWeight(iloc, Location(lat, lon, elev));
                  // The diagonal is 1, since the distance from a point to itself
                  // is 0, therefore its weight is 1.
                  matrix[ii][ii] = 1;
                  // The matrix is symmetric, so only compute one of the halves
                  for(int jj = ii+1; jj < N; jj++) {
                     Location jloc = validLocations[jj];
                     // Improve conditioning of matrix when you have two or more stations
                     // that are very close
                     float factor = 0.414 / 0.5;
                     float R = calcWeight(iloc, jloc)*factor;
                     // Store the number in both halves
                     matrix[ii][jj] = R;
                     matrix[jj][ii] = R;
                  }
               }

               // Compute (matrix)^-1
               vec2 inverse = Util::inverse(matrix);

               // Compute weights (matrix-vector product)
               std::vector<float> weights;
               weights.resize(N, 0);
               for(int ii = 0; ii < N; ii++) {
                  for(int jj = 0; jj < N; jj++) {
                     weights[ii] += inverse[ii][jj] * S[ii];
                  }
               }

               // Compute final bias (dot product of bias and weights)
               float finalBias = 0;
               for(int ii = 0; ii < N; ii++) {
                  float currBias = bias[ii];
                  if(!Util::isValid(currBias)) {
                     finalBias = Util::MV;
                     break;
                  }
                  finalBias += bias[ii]*weights[ii];
               }

               if(Util::isValid(finalBias)) {
                  // Apply bias to each ensemble member
                  for(int e = 0; e < nEns; e++) {
                     float rawValue = (*field)(i,j,e);
                     // Determine if we should turn kriging off based on auxillary variable
                     bool turnOn = Util::isValid(rawValue);
                     if(mAuxVariable != Variable::None) {
                        float aux = (*auxField)(i,j,e);
                        if(Util::isValid(aux)) {
                           if(aux < mLowerThreshold || aux > mUpperThreshold) {
                              turnOn = false;
                           }
                        }
                     }

                     if(turnOn) {
                        if(mOperator == Util::OperatorAdd) {
                           (*field)(i,j,e) += finalBias;
                        }
                        else if(mOperator == Util::OperatorSubtract) {
                           (*field)(i,j,e) -= finalBias;
                        }
                        else if(mOperator == Util::OperatorMultiply) {
                           // TODO: How do we ensure that the matrix is positive definite in this
                           // case?
                           (*field)(i,j,e) *= (finalBias+1);
                        }
                        else if(mOperator == Util::OperatorDivide) {
                           // TODO: How do we ensure that the matrix is positive definite in this
                           // case?
                           (*field)(i,j,e) /= (finalBias+1);
                        }
                        else {
                           Util::error("Unrecognized operator in CalibratorKriging");
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return true;
}

float CalibratorKriging::calcWeight(const Location& loc1, const Location& loc2) const {
   float weight = Util::MV;
   if(Util::isValid(loc1.lat()) && Util::isValid(loc1.lon()) && Util::isValid(loc2.lat()) && Util::isValid(loc2.lon()) && Util::isValid(loc1.elev()) && Util::isValid(loc2.elev())) {
      float horizDist = Util::getDistance(loc1.lat(), loc1.lon(), loc2.lat(), loc2.lon(), mUseApproxDistance);
      float vertDist  = fabs(loc1.elev() - loc2.elev());
      if(horizDist >= mRadius || vertDist >= mMaxElevDiff) {
         return 0;
      }
      if(mKrigingType == TypeCressman) {
         if(horizDist > mEfoldDist || vertDist > mMaxElevDiff)
            return 0;
         float horizWeight = (mEfoldDist*mEfoldDist - horizDist*horizDist) / (mEfoldDist*mEfoldDist + horizDist*horizDist);
         float vertWeight  = (mMaxElevDiff*mMaxElevDiff - vertDist*vertDist) / (mMaxElevDiff*mMaxElevDiff + vertDist*vertDist);
         weight = horizWeight * vertWeight;
      }
      else if(mKrigingType == TypeBarnes) {
         float horizWeight = exp(-horizDist*horizDist/(2*mEfoldDist*mEfoldDist));
         float vertWeight  = exp(-vertDist*vertDist/(2*mMaxElevDiff*mMaxElevDiff));
         weight = horizWeight * vertWeight;
         return weight;
      }
   }
   return weight;
}

Parameters CalibratorKriging::train(const TrainingData& iData, int iOffset) const {
   double timeStart = Util::clock();
   std::vector<ObsEns> data = iData.getData(iOffset);
   float totalBias = 0;
   int counter = 0;
   // Compute predictors in model
   for(int i = 0; i < data.size(); i++) {
      float obs = data[i].first;
      std::vector<float> ens = data[i].second;
      float mean = Util::calculateStat(ens, Util::StatTypeMean);
      if(Util::isValid(obs) && Util::isValid(mean)) {
         totalBias += obs - mean;
         counter++;
      }
   }

   float bias = 0;
   if(counter <= 0) {
      std::stringstream ss;
      ss << "CalibratorKriging: No valid data, no correction will be made.";
      Util::warning(ss.str());
      bias = 0;
   }
   else {
      bias = totalBias / counter;
   }

   std::vector<float> values;
   values.push_back(bias);
   Parameters par(values);

   double timeEnd = Util::clock();
   std::cout << "Time: " << timeEnd - timeStart << std::endl;
   return par;
}

std::string CalibratorKriging::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-c kriging","Spreads bias in space by using optimal interpolation. A parameter file is required, which must have one column with the bias.")<< std::endl;
   ss << Util::formatDescription("   radius=30000","Only use values from locations within this radius (in meters). Must be >= 0.") << std::endl;
   ss << Util::formatDescription("   efoldDist=30000","How fast should the weight of a station reduce with distance? For cressman: linearly decrease to this distance (in meters); For barnes: reduce to 1/e after this distance (in meters). Must be >= 0.") << std::endl;
   ss << Util::formatDescription("   maxElevDiff=100","What is the maximum elevation difference (in meters) that bias can be spread to? Must be >= 0.") << std::endl;
   ss << Util::formatDescription("   auxVariable=undef","Should an auxilary variable be used to turn off kriging? For example turn off kriging where there is precipitation.") << std::endl;
   ss << Util::formatDescription("   range=undef","What range of the auxillary variable should kriging be turned on for? For example use 0,0.3 to turn kriging off for precip > 0.3.") << std::endl;
   ss << Util::formatDescription("   type=cressman","Weighting function used in kriging. One of 'cressman', or 'barnes'.") << std::endl;
   ss << Util::formatDescription("   operator=add","How should the bias be applied to the raw forecast? One of 'add', 'subtract', 'multiply', 'divide'.") << std::endl;
   return ss.str();
}
