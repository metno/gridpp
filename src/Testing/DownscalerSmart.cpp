#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

namespace {
   class TestDownscalerSmart : public ::testing::Test {
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
         // iElev must be a flat array with lon varying fastest
         void setLatLonElev(FileFake& iFile, float iLat[], float iLon[], float iElev[]) {
            vec2 lat;
            vec2 lon;
            vec2 elev;
            int nLat = iFile.getNumLat(); 
            int nLon = iFile.getNumLon();
            lat.resize(nLat);
            lon.resize(nLat);
            elev.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               lat[i].resize(nLon);
               lon[i].resize(nLon);
               elev[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  lat[i][j] = iLat[i];
                  lon[i][j] = iLon[j];
                  elev[i][j] = iElev[j+i*nLon];
               }
            }
            iFile.setLats(lat);
            iFile.setLons(lon);
            iFile.setElevs(elev);
         };
      protected:
   };

   TEST_F(TestDownscalerSmart, isValid) {
      FileFake from(3,2,1,1);
      FileFake to(1,1,1,1);
      setLatLonElev(from, (float[]) {50,55,60}, (float[]){0,10}, (float[]){3, 15, 6, 30, 20, 11});
      setLatLonElev(to,   (float[]) {54},   (float[]){9}, (float[]){10});
      int searchRadius = 10; // Search the whole grid
      int numSmart = 2;

      vec3Int I, J;
      DownscalerSmart::getSmartNeighbours(from, to, searchRadius, numSmart, I, J);

      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, I[0].size());
      ASSERT_EQ(2, I[0][0].size());

      EXPECT_EQ(2, I[0][0][0]);
      EXPECT_EQ(1, J[0][0][0]);
      EXPECT_EQ(1, I[0][0][1]);
      EXPECT_EQ(0, J[0][0][1]);
   }
   TEST_F(TestDownscalerSmart, 10x10) {
      DownscalerSmart d(Variable::T);
      d.setSearchRadius(3);
      d.setNumSmart(2);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      FileFake to(1,4,1,1);
      // Case 1: West boundary outside domain
      // Case 2: Within domain
      // Case 3/4: Nearest neighbour is on the boundary, so only the western half is used
      setLatLonElev(to, (float[]) {5}, (float[]){2,3,12,20}, (float[]){120, 80, 600, 600});
      d.downscale(from, to);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT.size());
      ASSERT_EQ(4, toT[0].size());
      EXPECT_FLOAT_EQ(303,   toT[0][0][0]);
      EXPECT_FLOAT_EQ(304.5, toT[0][1][0]);
      EXPECT_FLOAT_EQ(305.5, toT[0][2][0]);
      EXPECT_FLOAT_EQ(305.5,   toT[0][3][0]);
      vec3Int I, J;
      d.getSmartNeighbours(from, to, I, J);
      EXPECT_EQ(4, I[0][2][0]);
      EXPECT_EQ(9, J[0][2][0]);
      EXPECT_EQ(3, I[0][2][1]);
      EXPECT_EQ(8, J[0][2][1]);
   }
   TEST_F(TestDownscalerSmart, fewerPointsThanSmart) {
      FileFake from(5,3,1,1);
      setLatLonElev(from, (float[]) {4,5,6,7,8}, (float[]){5,10,15}, (float[]){70,50,20,80,60,70,50,40,30,20,10,40,50,30,60});
      FileFake to(1,4,1,1);
      setLatLonElev(to, (float[]) {5.5}, (float[]){2,4, 10,20}, (float[]){120, 80, 600, 600});
      vec3Int I, J;
      DownscalerSmart::getSmartNeighbours(from, to, 1, 20, I, J);
      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, J.size());
      ASSERT_EQ(4, I[0].size());
      ASSERT_EQ(4, J[0].size());
      EXPECT_EQ(6, I[0][0].size());
      EXPECT_EQ(6, J[0][0].size());
      EXPECT_EQ(6, I[0][1].size());
      EXPECT_EQ(6, J[0][1].size());
      EXPECT_EQ(9, I[0][2].size());
      EXPECT_EQ(9, J[0][2].size());
      EXPECT_EQ(6, I[0][3].size());
      EXPECT_EQ(6, J[0][3].size());
   }
   TEST_F(TestDownscalerSmart, setGetSearchRadius) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      DownscalerSmart d(Variable::Precip);
      // Check that default is valid
      EXPECT_GE(d.getSearchRadius(), 0);
      d.setSearchRadius(5);
      EXPECT_FLOAT_EQ(5, d.getSearchRadius());
      d.setSearchRadius(0);
      EXPECT_FLOAT_EQ(0, d.getSearchRadius());

      // Invalid values
      EXPECT_DEATH(d.setSearchRadius(-1), ".*");
      EXPECT_DEATH(d.setSearchRadius(Util::MV), ".*");
   }
   TEST_F(TestDownscalerSmart, setGetNumSmart) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      DownscalerSmart d(Variable::Precip);
      // Check that default is valid
      EXPECT_GT(d.getNumSmart(), 0);
      d.setNumSmart(5);
      EXPECT_FLOAT_EQ(5, d.getNumSmart());
      d.setNumSmart(1);
      EXPECT_FLOAT_EQ(1, d.getNumSmart());

      // Invalid values
      EXPECT_DEATH(d.setNumSmart(-1), ".*");
      EXPECT_DEATH(d.setNumSmart(0), ".*");
      EXPECT_DEATH(d.setNumSmart(Util::MV), ".*");
   }
   TEST_F(TestDownscalerSmart, getNumSearchPoints) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      DownscalerSmart d(Variable::Precip);
      // Check that default is valid
      EXPECT_GT(d.getNumSearchPoints(), 0);
      d.setNumSmart(5);
      d.setSearchRadius(5);
      EXPECT_FLOAT_EQ(121, d.getNumSearchPoints());
      EXPECT_FLOAT_EQ(121, DownscalerSmart::getNumSearchPoints(5));
      d.setNumSmart(1);
      EXPECT_FLOAT_EQ(121, d.getNumSearchPoints());
      d.setSearchRadius(3);
      EXPECT_FLOAT_EQ(49, d.getNumSearchPoints());
      EXPECT_FLOAT_EQ(49, DownscalerSmart::getNumSearchPoints(3));
   }
   TEST_F(TestDownscalerSmart, description) {
      DownscalerSmart::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
