#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class TestDownscalerNearestNeighbour : public ::testing::Test {
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
         void setLatLon(FileFake& iFile, float iLat[], float iLon[]) {
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
   };

   TEST_F(TestDownscalerNearestNeighbour, description) {
      DownscalerNearestNeighbour::description();
   }
   TEST_F(TestDownscalerNearestNeighbour, downscale) {
      DownscalerNearestNeighbour d(Variable::T);
      FileFake from(3,2,1,1);
      FileFake to(2,2,1,1);
      setLatLon(from, (float[]) {60,50,55}, (float[]){5,4});
      setLatLon(to,   (float[]) {56,49},    (float[]){3,4.6});
      d.downscale(from, to);
      const Field& fromT = *from.getField(Variable::T, 0);
      const Field& toT   = *to.getField(Variable::T, 0);
      EXPECT_FLOAT_EQ(fromT(2,1,0), toT(0,0,0));
      EXPECT_FLOAT_EQ(fromT(2,0,0), toT(0,1,0));
      EXPECT_FLOAT_EQ(fromT(1,1,0), toT(1,0,0));
      EXPECT_FLOAT_EQ(fromT(1,0,0), toT(1,1,0));
   }
   TEST_F(TestDownscalerNearestNeighbour, 10x10) {
      DownscalerNearestNeighbour d(Variable::T);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      FileFake to(2,2,1,from.getNumTime());
      setLatLon(to,   (float[]) {5,11},    (float[]){2,3});
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(2, toT.getNumLat());
      ASSERT_EQ(2, toT.getNumLon());
      EXPECT_FLOAT_EQ(301, toT(0,0,0));
      EXPECT_FLOAT_EQ(302, toT(0,1,0));
      EXPECT_FLOAT_EQ(309, toT(1,0,0));
      EXPECT_FLOAT_EQ(301, toT(1,1,0));
      vec2Int I, J;
      d.getNearestNeighbour(from, to, I, J);
      EXPECT_EQ(9, I[1][1]);
      EXPECT_EQ(3, J[1][1]);
      EXPECT_FLOAT_EQ(301, fromT(9,3,0));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
