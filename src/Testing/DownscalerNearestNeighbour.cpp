#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

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

   TEST_F(TestDownscalerNearestNeighbour, isValid) {
      FileFake from(3,2,1,1);
      FileFake to(2,2,1,1);
      setLatLon(from, (float[]) {50,55,60}, (float[]){0,10});
      setLatLon(to,   (float[]) {40, 54},   (float[]){-1,9});

      vec2Int I, J, If, Jf;
      DownscalerNearestNeighbour::getNearestNeighbour(from, to, I, J);
      DownscalerNearestNeighbour::getNearestNeighbourFast(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());

      EXPECT_EQ(0, I[0][0]);
      EXPECT_EQ(0, J[0][0]);
      EXPECT_EQ(0, I[0][1]);
      EXPECT_EQ(1, J[0][1]);
      EXPECT_EQ(1, I[1][0]);
      EXPECT_EQ(0, J[1][0]);
      EXPECT_EQ(1, I[1][1]);
      EXPECT_EQ(1, J[1][1]);
   }
   TEST_F(TestDownscalerNearestNeighbour, emptyFrom) {
      FileFake from(0,0,1,1);
      FileFake to(1,1,1,1);
      setLatLon(to,   (float[]) {40},   (float[]){9});

      vec2Int I, J, If, Jf;
      DownscalerNearestNeighbour::getNearestNeighbour(from, to, I, J);
      DownscalerNearestNeighbour::getNearestNeighbourFast(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, I[0].size());

      EXPECT_EQ(Util::MV, I[0][0]);
      EXPECT_EQ(Util::MV, J[0][0]);
   }
   TEST_F(TestDownscalerNearestNeighbour, emptyTo) {
      FileFake from(3,2,1,1);
      FileFake to(0,0,1,1);
      setLatLon(from, (float[]) {50,55,60}, (float[]){0,10});

      vec2Int I, J, If, Jf;
      DownscalerNearestNeighbour::getNearestNeighbour(from, to, I, J);
      DownscalerNearestNeighbour::getNearestNeighbourFast(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(0, I.size());
   }
   TEST_F(TestDownscalerNearestNeighbour, missingLatLon) {
      FileFake from(3,2,1,1);
      FileFake to(2,2,1,1);
      setLatLon(from, (float[]) {50,Util::MV,60}, (float[]){0,Util::MV});
      setLatLon(to,   (float[]) {40, 54},   (float[]){-1,Util::MV});

      vec2Int I, J, If, Jf;
      DownscalerNearestNeighbour::getNearestNeighbour(from, to, I, J);
      DownscalerNearestNeighbour::getNearestNeighbourFast(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());

      EXPECT_EQ(0, I[0][0]);
      EXPECT_EQ(0, J[0][0]);
      EXPECT_EQ(Util::MV, I[0][1]);
      EXPECT_EQ(Util::MV, J[0][1]);
      EXPECT_EQ(0, I[1][0]);
      EXPECT_EQ(0, J[1][0]);
      EXPECT_EQ(Util::MV, I[1][1]);
      EXPECT_EQ(Util::MV, J[1][1]);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
