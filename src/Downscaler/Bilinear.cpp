#include "Bilinear.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

// std::map<const File*, std::map<const File*, std::pair<vec2Int, vec2Int> > > DownscalerBilinear::mNeighbourCache;

DownscalerBilinear::DownscalerBilinear(Variable::Type iVariable, const Options& iOptions) :
      Downscaler(iVariable, iOptions) {
}

void DownscalerBilinear::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumLat();
   int nLon = iOutput.getNumLon();
   int nEns = iOutput.getNumEns();
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

      // Algorithm from here:
      // https://stackoverflow.com/questions/23920976/bilinear-interpolation-with-non-aligned-input-points
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int I = nearestI[i][j];
            int J = nearestJ[i][j];
            int I1 = I-1;
            int I2 = I;
            int J1 = J-1;
            int J2 = J;
            float lat = olats[i][j];
            float lon = olons[i][j];
            std::cout << ilons[I][J] << " " << ilons[I][J-1] << " " << ilons[I][J+1] << std::endl;
            std::cout << ilats[I][J] << " " << ilats[I][J-1] << " " << ilats[I][J+1] << std::endl;
#if 1
            if(ilons[I][J] < lon) {
               I1 = I;
               I2 = I + 1;
               std::cout << "Switch lon" << std::endl;
            }
            if(ilats[I][J] < lat) {
               J1 = J;
               J2 = J + 1;
               std::cout << "Switch lat" << std::endl;
            }
#endif
            float x1 = ilons[I1][J1];
            float x2 = ilons[I2][J1];
            float x3 = ilons[I1][J2];
            float x4 = ilons[I2][J2];
            float y1 = ilats[I1][J1];
            float y2 = ilats[I2][J1];
            float y3 = ilats[I1][J2];
            float y4 = ilats[I2][J2];
            std::cout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << lon << std::endl;
            std::cout << y1 << " " << y2 << " " << y3 << " " << y4 << " " << lat << std::endl;
            float a = -x1 + x3;
            float b = -x1 + x2;
            float c = x1 - x2 - x3 + x4;
            float d = lon - x1;
            float e = -y1 + y3;
            float f = -y1 + y2;
            float g = y1 - y2 - y3 + y4;
            float h = lat - y1;
            float alpha = -(b * e - a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
            float beta  = (b * e - a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
#if 0
            if(alpha < 0 || alpha > 1 || beta < 0 || beta > 1) {
               // alpha = (-b * e + a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
               // beta  = -((-b * e + a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g));
               alpha = -(b * e - a * f + d * g - c * h - sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
               beta  = (b * e - a * f - d * g + c * h - sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
            }
#endif
            std::cout << alpha << " " << beta << std::endl;

            for(int ee = 0; ee < nEns; ee++) {
               if(Util::isValid(I) && Util::isValid(J)) {
                  float v1 = ifield(I1, J1, ee);
                  float v2 = ifield(I2, J1, ee);
                  float v3 = ifield(I1, J2, ee);
                  float v4 = ifield(I2, J2, ee);
                  float value = (1 - alpha) * ((1 - beta) * v1 + beta * v2) + alpha * ((1 - beta) * v3 + beta * v4);
                  std::cout << v1 << " " << v2 << " " << v3 << " " << v4 << " " << value << std::endl;
                  ofield(i,j,ee) = value;
               }
               else
                  ofield(i,j,ee) = Util::MV;
            }
         }
      }
   }
}

float DownscalerBilinear::bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3) const {
   float a = -x0 + x2;
   float b = -x0 + x1;
   float c = x0 - x1 - x2 + x3;
   float d = x - x0;
   float e = -y0 + y2;
   float f = -y0 + y1;
   float g = y0 - y1 - y2 + y3;
   float h = y - y0;
   float alpha = -(b * e - a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
   float beta  = (b * e - a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
#if 0
   if(alpha < 0 || alpha > 1 || beta < 0 || beta > 1) {
      // alpha = (-b * e + a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
      // beta  = -((-b * e + a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g));
      alpha = -(b * e - a * f + d * g - c * h - sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
      beta  = (b * e - a * f - d * g + c * h - sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
      }
#endif
   float value = (1 - alpha) * ((1 - beta) * v0 + beta * v1) + alpha * ((1 - beta) * v2 + beta * v3);
   return value;
}

std::string DownscalerBilinear::description() {
   std::stringstream ss;
   ss << Util::formatDescription("-d bilinear", "Uses bilinear interpolation") << std::endl;
   return ss.str();
}
