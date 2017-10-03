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
            // Use the four points surrounding the lookup point. Use the nearest neighbour
            // and then figure out what side of this point the other there points are.
            int I = nearestI[i][j];
            int J = nearestJ[i][j];
            int I1 = I-1;
            int I2 = I;
            int J1 = J-1;
            int J2 = J;
            float lat = olats[i][j];
            float lon = olons[i][j];
            // std::cout << ilons[I][J] << " " << ilons[I][J-1] << " " << ilons[I][J+1] << std::endl;
            // std::cout << ilats[I][J] << " " << ilats[I][J-1] << " " << ilats[I][J+1] << std::endl;
            if(ilons[I][J] < lon) {
               I1 = I;
               I2 = I + 1;
            }
            if(ilats[I][J] < lat) {
               J1 = J;
               J2 = J + 1;
            }
            float x0 = ilons[I1][J1];
            float x1 = ilons[I2][J1];
            float x2 = ilons[I1][J2];
            float x3 = ilons[I2][J2];
            float y0 = ilats[I1][J1];
            float y1 = ilats[I2][J1];
            float y2 = ilats[I1][J2];
            float y3 = ilats[I2][J2];
            // std::cout << x0 << " " << x1 << " " << x2 << " " << x3 << " " << lon << std::endl;
            // std::cout << y0 << " " << y1 << " " << y2 << " " << y3 << " " << lat << std::endl;

            for(int ee = 0; ee < nEns; ee++) {
               if(Util::isValid(I) && Util::isValid(J)) {
                  float v0 = ifield(I1, J1, ee);
                  float v1 = ifield(I2, J1, ee);
                  float v2 = ifield(I1, J2, ee);
                  float v3 = ifield(I2, J2, ee);
                  float value = bilinear(lon, lat, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3);
                  // std::cout << v0 << " " << v1 << " " << v2 << " " << v3 << " " << value << std::endl;
                  ofield(i,j,ee) = value;
               }
               else
                  ofield(i,j,ee) = Util::MV;
            }
         }
      }
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
   // std::cout << y0 << " " << y1 << " " << y2 << " " << y3 << std::endl;
   // std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << std::endl;
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
