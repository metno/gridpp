#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class TestDownscalerBilinear : public ::testing::Test {
      public:
         vec2 makeVec2(int nLat, int nLon, const std::vector<float>& values) {
            vec2 grid;
            grid.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               grid[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  int index = i*nLon + j;
                  grid[i][j] = values[index];
               }
            }
            return grid;
         };
         void setLatLon(FileFake& iFile, const float iLat[], const float iLon[]) {
            vec2 lat;
            vec2 lon;
            int nLat = iFile.getNumLat();
            int nLon = iFile.getNumLon();
            lat.resize(nLat);
            lon.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               lat[i].resize(nLon);
               lon[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  lat[i][j] = iLat[i];
                  lon[i][j] = iLon[j];
               }
            }
            iFile.setLats(lat);
            iFile.setLons(lon);
         };
      protected:
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(TestDownscalerBilinear, description) {
      DownscalerBilinear::description();
   }
   TEST_F(TestDownscalerBilinear, square) {
      // Verifying against:
      // https://www.ajdesigner.com/phpinterpolation/bilinear_interpolation_equation.php
      // Using values:
      // 1.1, 1.3, 2.2, 4.8, 0.6, 0.3, 10.3, 10.9, 10, 10.7
      float x0 = 1.1;
      float x1 = 2.2;
      float y0 = 0.3;
      float y1 = 4.8;
      float v0 = 10;
      float v1 = 10.3;
      float v2 = 10.7;
      float v3 = 10.9;
      float value = DownscalerBilinear::bilinear(1.3, 0.6,
            x0, x0, x1, x1,
            y0, y1, y0, y1,
            v0, v1, v2, v3);
      EXPECT_FLOAT_EQ(10.146060606061, value);
      // Test corners
      EXPECT_FLOAT_EQ(10, DownscalerBilinear::bilinear(1.1, 0.3, 1.1, 1.1, 2.2, 2.2, 0.3, 4.8, 0.3, 4.8, 10,10.3,10.7,10.9));
      EXPECT_FLOAT_EQ(10.3, DownscalerBilinear::bilinear(1.1, 4.8, 1.1, 1.1, 2.2, 2.2, 0.3, 4.8, 0.3, 4.8, 10,10.3,10.7,10.9));
      EXPECT_FLOAT_EQ(10.7, DownscalerBilinear::bilinear(2.2, 0.3, 1.1, 1.1, 2.2, 2.2, 0.3, 4.8, 0.3, 4.8, 10,10.3,10.7,10.9));
      EXPECT_FLOAT_EQ(10.9, DownscalerBilinear::bilinear(2.2, 4.8, 1.1, 1.1, 2.2, 2.2, 0.3, 4.8, 0.3, 4.8, 10,10.3,10.7,10.9));

      // Now reverse the points, so that v2,v3 are to the left of v0,v1
      EXPECT_FLOAT_EQ(10.146060606061, DownscalerBilinear::bilinear(1.3,0.6,x1,x1,x0,x0,y0,y1,y0,y1,v2,v3,v0,v1));
      // Put v1,v3 below v0,v1
      EXPECT_FLOAT_EQ(10.146060606061, DownscalerBilinear::bilinear(1.3,0.6,x0,x0,x1,x1,y1,y0,y1,y0,v1,v0,v3,v2));
   }
   TEST_F(TestDownscalerBilinear, downscale) {
      // Same as above, but using an input file
      DownscalerBilinear d(mVariable, mVariable, Options());
      FileFake from(Options("nLat=2 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {0.3, 4.8}, (const float[]){1.1, 2.2});
      setLatLon(to,   (const float[]) {0.6},    (const float[]){1.3});

      Field& fromT = *from.getField(mVariable, 0);
      fromT(0,0,0) = 10;
      fromT(1,0,0) = 10.3;
      fromT(0,1,0) = 10.7;
      fromT(1,1,0) = 10.9;
      d.downscale(from, to);
      const Field& toT   = *to.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(10.14606060606061, toT(0,0,0));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
