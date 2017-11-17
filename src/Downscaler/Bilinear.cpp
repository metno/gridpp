#include "Bilinear.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerBilinear::DownscalerBilinear(Variable::Type iVariable, const Options& iOptions) :
      Downscaler(iVariable, iOptions) {
}

void DownscalerBilinear::downscaleCore(const File& iInput, File& iOutput) const {
   int nTime = iInput.getNumTime();
   vec2 ilats = iInput.getLats();
   vec2 ilons = iInput.getLons();
   vec2 olats = iOutput.getLats();
   vec2 olons = iOutput.getLons();

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   Downscaler::getNearestNeighbour(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mVariable, t);
      Field& ofield = *iOutput.getField(mVariable, t);
      downscaleField(ifield, ofield, ilats, ilons, olats, olons, nearestI, nearestJ);
   }
}

float DownscalerBilinear::bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3) {
   float a = -x0 + x2;
   float b = -x0 + x1;
   float c = x0 - x1 - x2 + x3;
   float d = x - x0;
   float e = -y0 + y2;
   float f = -y0 + y1;
   float g = y0 - y1 - y2 + y3;
   float h = y - y0;
   float alpha = 0;
   float beta  = 0;
   if(c*e == a*g || c*f == b*g) {
      alpha = (x - x0) / (x2 - x0);
      beta = (y - y0) / (y1 - y0);
      assert(x0 == x1);
      assert(x2 == x3);
      assert(y0 == y2);
      assert(y1 == y3);
   }
   else {
      assert(2 * c * e - 2 * a * g != 0);
      assert(2 * c * f - 2 * b * g != 0);
      alpha = -(b * e - a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
      beta  = (b * e - a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
   }
   float value = (1 - alpha) * ((1 - beta) * v0 + beta * v1) + alpha * ((1 - beta) * v2 + beta * v3);
   return value;
}

vec2 DownscalerBilinear::downscaleVec(const vec2& iInput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ) {

   int nLat = iOutputLats.size();
   int nLon = iOutputLats[0].size();

   vec2 output;
   output.resize(nLat);
   for(int i = 0; i < nLat; i++)
      output[i].resize(nLon);

   // Algorithm from here:
   // https://stackoverflow.com/questions/23920976/bilinear-interpolation-with-non-aligned-input-points
   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         // Use the four points surrounding the lookup point. Use the nearest neighbour
         // and then figure out what side of this point the other there points are.
         int I = nearestI[i][j];
         int J = nearestJ[i][j];
         float lat = iOutputLats[i][j];
         float lon = iOutputLons[i][j];
         // Find out if current lookup point is right/above the nearest neighbour
         bool isRight = lon > iInputLons[I][J];
         bool isAbove = lat > iInputLats[I][J];
         int I1 = Util::MV;
         int I2 = Util::MV;
         int J1 = Util::MV;
         int J2 = Util::MV;
         // Find out which Is to use
         if((isAbove && iInputLats[I+1][J] > iInputLats[I][J]) || (!isAbove && iInputLats[I+1][J] < iInputLats[I][J])) {
            I1 = I;
            I2 = I + 1;
         }
         else {
            I1 = I - 1;
            I2 = I;
         }
         // Find out which Js to use
         if((isRight && iInputLons[I][J+1] > iInputLons[I][J]) || (!isRight && iInputLons[I][J+1] < iInputLons[I][J])) {
            J1 = J;
            J2 = J + 1;
         }
         else {
            J1 = J - 1;
            J2 = J;
         }
         if(I2 < iInputLons.size() && J2 < iInputLons[I].size()) {
            float x0 = iInputLons[I1][J1];
            float x1 = iInputLons[I2][J1];
            float x2 = iInputLons[I1][J2];
            float x3 = iInputLons[I2][J2];
            float y0 = iInputLats[I1][J1];
            float y1 = iInputLats[I2][J1];
            float y2 = iInputLats[I1][J2];
            float y3 = iInputLats[I2][J2];

            float v0 = iInput[I1][J1];
            float v1 = iInput[I2][J1];
            float v2 = iInput[I1][J2];
            float v3 = iInput[I2][J2];
            if(Util::isValid(v0) && Util::isValid(v1) &&Util::isValid(v2) &&Util::isValid(v3)) {
               float value = bilinear(lon, lat, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3);
               output[i][j] = value;
            }
            else
               output[i][j] = iInput[I][J];
         }
         else {
            // The point is outside the input domain. Revert to nearest neighbour
            // Util::warning("Point is outside domain, cannot bilinearly interpolate");
            output[i][j] = iInput[I][J];
         }
      }
   }
   return output;
}

void DownscalerBilinear::downscaleField(const Field& iInput, Field& iOutput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ) {

   int nLat = iOutput.getNumLat();
   int nLon = iOutput.getNumLon();
   int nEns = iOutput.getNumEns();

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         int I = nearestI[i][j];
         int J = nearestJ[i][j];
         float lat = iOutputLats[i][j];
         float lon = iOutputLons[i][j];
         // Find out if current lookup point is right/above the nearest neighbour
         bool isRight = lon > iInputLons[I][J];
         bool isAbove = lat > iInputLats[I][J];
         int I1 = Util::MV;
         int I2 = Util::MV;
         int J1 = Util::MV;
         int J2 = Util::MV;
         // Find out which Is to use
         if((isAbove && iInputLats[I+1][J] > iInputLats[I][J]) || (!isAbove && iInputLats[I+1][J] < iInputLats[I][J])) {
            I1 = I;
            I2 = I + 1;
         }
         else {
            I1 = I - 1;
            I2 = I;
         }
         // Find out which Js to use
         if((isRight && iInputLons[I][J+1] > iInputLons[I][J]) || (!isRight && iInputLons[I][J+1] < iInputLons[I][J])) {
            J1 = J;
            J2 = J + 1;
         }
         else {
            J1 = J - 1;
            J2 = J;
         }
         if(I2 < iInputLons.size() && J2 < iInputLons[I].size()) {
            float x0 = iInputLons[I1][J1];
            float x1 = iInputLons[I2][J1];
            float x2 = iInputLons[I1][J2];
            float x3 = iInputLons[I2][J2];
            float y0 = iInputLats[I1][J1];
            float y1 = iInputLats[I2][J1];
            float y2 = iInputLats[I1][J2];
            float y3 = iInputLats[I2][J2];

            for(int ee = 0; ee < nEns; ee++) {
               float v0 = iInput(I1, J1, ee);
               float v1 = iInput(I2, J1, ee);
               float v2 = iInput(I1, J2, ee);
               float v3 = iInput(I2, J2, ee);
               if(Util::isValid(v0) && Util::isValid(v1) &&Util::isValid(v2) &&Util::isValid(v3)) {
                  float value = bilinear(lon, lat, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3);
                  iOutput(i,j,ee) = value;
               }
               else
                  iOutput(i,j,ee) = Util::MV;
            }
         }
         else {
            // The point is outside the input domain. Revert to nearest neighbour
            // Util::warning("Point is outside domain, cannot bilinearly interpolate");
            for(int ee = 0; ee < nEns; ee++) {
               iOutput(i,j,ee) = iInput(I, J, ee);
            }
         }
      }
   }
}

std::string DownscalerBilinear::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-d bilinear", "Bilinear interpolation using the 4 points surrounding the lookup point. If the lookup point is outside the input domain, then the nearest neighbour is used.") << std::endl;
   return ss.str();
}
