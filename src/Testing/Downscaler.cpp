#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>
#include <boost/uuid/uuid.hpp>

namespace {
   class TestDownscaler : public ::testing::Test {
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

   TEST_F(TestDownscaler, validDownscalers) {
      Downscaler* d0 = Downscaler::getScheme("nearestNeighbour", Variable::T, Options());
      Downscaler* d1 = Downscaler::getScheme("smart", Variable::T, Options("searchRadius=3 numSmart=2 minElevDiff=400"));
      Downscaler* d2 = Downscaler::getScheme("gradient", Variable::T, Options("searchRadius=5 constantGradient=0.04"));
      EXPECT_EQ(3, ((DownscalerSmart*) d1)->getSearchRadius());
      EXPECT_EQ(2, ((DownscalerSmart*) d1)->getNumSmart());
      EXPECT_EQ(400, ((DownscalerSmart*) d1)->getMinElevDiff());
      EXPECT_EQ(5, ((DownscalerGradient*) d2)->getSearchRadius());
      EXPECT_FLOAT_EQ(0.04, ((DownscalerGradient*) d2)->getConstantGradient());
   }
   TEST_F(TestDownscaler, invalidDownscalers) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(Downscaler::getScheme("woehiowciwofew", Variable::T, Options("searchRadius=5")), ".*");
      EXPECT_DEATH(Downscaler::getScheme("woehiowciwofew", Variable::T, Options()), ".*");
   }
   TEST_F(TestDownscaler, validNearestNeighbours) {
      FileFake from(3,2,1,1);
      FileFake to(2,2,1,1);
      setLatLon(from, (float[]) {50,55,60}, (float[]){0,10});
      setLatLon(to,   (float[]) {40, 54.99},   (float[]){-1,9.99});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::getNearestNeighbourFast(from, to, If, Jf);
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
   TEST_F(TestDownscaler, missingLatLon) {
      FileFake from(3,2,1,1);
      FileFake to(2,2,1,1);
      setLatLon(from, (float[]) {50,Util::MV,60}, (float[]){0,Util::MV});
      setLatLon(to,   (float[]) {40, 54},   (float[]){-1,Util::MV});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::getNearestNeighbourFast(from, to, If, Jf);
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
   TEST_F(TestDownscaler, identicalGrid) {
      int nLat = 3;
      int nLon = 4;
      FileFake from(nLat,nLon,1,1);
      FileFake to(nLat,nLon,1,1);
      setLatLon(from, (float[]) {50,55,60}, (float[]){0,3,5,10});
      setLatLon(to,   (float[]) {50,55,60}, (float[]){0,3,5,10});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::getNearestNeighbourFast(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(nLat, I.size());
      ASSERT_EQ(nLon, I[0].size());

      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            EXPECT_EQ(i, I[i][j]);
            EXPECT_EQ(j, J[i][j]);
         }
      }
   }
   TEST_F(TestDownscaler, unsortedGrid) {
      FileFake from(3,2,1,1);
      FileFake to(2,2,1,1);
      setLatLon(from, (float[]) {60,50,55}, (float[]){5,4});
      setLatLon(to,   (float[]) {56,49},    (float[]){3,4.6});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::getNearestNeighbour(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());

      EXPECT_EQ(2, I[0][0]);
      EXPECT_EQ(1, J[0][0]);
      EXPECT_EQ(2, I[0][0]);
      EXPECT_EQ(0, J[0][1]);
      EXPECT_EQ(1, I[1][0]);
      EXPECT_EQ(1, J[1][0]);
      EXPECT_EQ(1, I[1][0]);
      EXPECT_EQ(0, J[1][1]);
   }
   TEST_F(TestDownscaler, cache) {
      FileFake from(3,2,1,1);
      FileFake to(2,2,1,1);
      setLatLon(from, (float[]) {60,50,55}, (float[]){5,4});
      setLatLon(to,   (float[]) {56,49},    (float[]){3,4.6});

      vec2Int I, J;
      Downscaler::getNearestNeighbour(from, to, I, J);
      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());
      EXPECT_EQ(2, I[0][0]);
      EXPECT_EQ(1, J[0][0]);

      // Don't change the grid
      boost::uuids::uuid idBefore = from.getUniqueTag();
      setLatLon(from, (float[]) {60,50,55}, (float[]){5,4});
      boost::uuids::uuid idAfter = from.getUniqueTag();
      EXPECT_EQ(idBefore, idAfter);
      Downscaler::getNearestNeighbour(from, to, I, J);
      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());
      EXPECT_EQ(2, I[0][0]);
      EXPECT_EQ(1, J[0][0]);

      // Change the grid
      idBefore = from.getUniqueTag();
      setLatLon(from, (float[]) {60,55,50}, (float[]){5,4});
      idAfter = from.getUniqueTag();
      EXPECT_NE(idBefore, idAfter);
      Downscaler::getNearestNeighbour(from, to, I, J);
      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());
      EXPECT_EQ(1, I[0][0]);
      EXPECT_EQ(1, J[0][0]);
   }
   TEST_F(TestDownscaler, copyConstructor) {
      FileFake from(3,2,1,1);
      FileFake to = from;
      vec2Int I, J;
      EXPECT_EQ(from.getUniqueTag(), to.getUniqueTag());
      Downscaler::getNearestNeighbour(from, to, I, J);
      ASSERT_EQ(3, I.size());
      ASSERT_EQ(2, I[0].size());
      for(int i = 0; i < I.size(); i++) {
         for(int j = 0; j < I[0].size(); j++) {
            EXPECT_EQ(i, I[i][j]);
            EXPECT_EQ(j, J[i][j]);
         }
      }
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
