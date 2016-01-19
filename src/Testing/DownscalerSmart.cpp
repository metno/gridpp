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
         void setLatLonElev(FileFake& iFile, const float iLat[], const float iLon[], const float iElev[]) {
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
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      setLatLonElev(from, (const float[]) {50,55,60}, (const float[]){0,10}, (const float[]){3, 15, 6, 30, 20, 11});
      setLatLonElev(to,   (const float[]) {54},   (const float[]){9}, (const float[]){10});

      vec3Int I, J;
      DownscalerSmart d(Variable::Precip, Options());
      d.setSearchRadius(10); // Search the whole grid
      d.setNumSmart(2);
      d.getSmartNeighbours(from, to, I, J);

      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, I[0].size());
      ASSERT_EQ(2, I[0][0].size());

      EXPECT_EQ(2, I[0][0][0]);
      EXPECT_EQ(1, J[0][0][0]);
      EXPECT_EQ(1, I[0][0][1]);
      EXPECT_EQ(0, J[0][0][1]);
   }
   TEST_F(TestDownscalerSmart, 10x10) {
      DownscalerSmart d(Variable::T, Options());
      d.setSearchRadius(3);
      d.setNumSmart(2);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      std::stringstream ss;
      ss << "nLat=1 nLon=5 nEns=1 nTime=" << from.getNumTime();
      FileFake to(Options(ss.str()));
      // Case 1: West boundary outside domain
      // Case 2: Within domain
      // Case 3/4: Nearest neighbour is on the boundary, so only the western half is used
      // Case 5: Elevation is missing, so use the nearest neighbour
      float elev[] = {120, 80, 600, 600, Util::MV};
      setLatLonElev(to, (const float[]) {5}, (const float[]){2,3,12,20,2}, elev);
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT.getNumLat());
      ASSERT_EQ(5, toT.getNumLon());
      EXPECT_FLOAT_EQ(303,   toT(0,0,0));
      EXPECT_FLOAT_EQ(304.5, toT(0,1,0));
      EXPECT_FLOAT_EQ(305.5, toT(0,2,0));
      EXPECT_FLOAT_EQ(305.5, toT(0,3,0));
      EXPECT_FLOAT_EQ(301,   toT(0,4,0));
      vec3Int I, J;
      d.getSmartNeighbours(from, to, I, J);
      EXPECT_EQ(4, I[0][2][0]);
      EXPECT_EQ(9, J[0][2][0]);
      EXPECT_EQ(3, I[0][2][1]);
      EXPECT_EQ(8, J[0][2][1]);

      ASSERT_EQ(1, I[0][4].size());
      ASSERT_EQ(1, J[0][4].size());
      EXPECT_EQ(5, I[0][4][0]);
      EXPECT_EQ(2, J[0][4][0]);
   }
   TEST_F(TestDownscalerSmart, 10x10minElevDiff) {
      DownscalerSmart d(Variable::T, Options());
      d.setSearchRadius(3);
      d.setNumSmart(2);
      d.setMinElevDiff(109);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      std::stringstream ss;
      ss << "nLat=1 nLon=3 nEns=1 nTime=" << from.getNumTime();
      FileFake to(Options(ss.str()));
      float elev[] = {120, 50, Util::MV};
      setLatLonElev(to, (const float[]) {5}, (const float[]){2,2,2}, elev);
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT.getNumLat());
      ASSERT_EQ(3, toT.getNumLon());
      EXPECT_FLOAT_EQ(301, toT(0,0,0)); // Nearest neighbour
      EXPECT_FLOAT_EQ(304, toT(0,1,0));
      EXPECT_FLOAT_EQ(301, toT(0,2,0)); // Nearest neighbour
      vec3Int I, J;
      d.getSmartNeighbours(from, to, I, J);
      ASSERT_EQ(2, I[0][1].size());
      ASSERT_EQ(2, J[0][1].size());
      EXPECT_EQ(4, I[0][1][0]);
      EXPECT_EQ(5, J[0][1][0]);
      EXPECT_EQ(2, I[0][1][1]);
      EXPECT_EQ(4, J[0][1][1]);

      // All nearest neighbours
      d.setMinElevDiff(1200);
      status = d.downscale(from, to);
      EXPECT_TRUE(status);
      EXPECT_FLOAT_EQ(301, toT(0,0,0));
      EXPECT_FLOAT_EQ(301, toT(0,1,0));
      EXPECT_FLOAT_EQ(301, toT(0,2,0));
   }
   TEST_F(TestDownscalerSmart, minElevDiff) {
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      setLatLonElev(from, (const float[]) {50,55,60}, (const float[]){0,10}, (const float[]){3, 15, 6, 30, 20, 11});
      setLatLonElev(to,   (const float[]) {55},       (const float[]){10},   (const float[]){22.3});

      DownscalerSmart d(Variable::T, Options());
      d.setSearchRadius(3);
      d.setNumSmart(2);
      vec3Int I, J;

      // Use nearest neighbour (at 30m)
      d.setMinElevDiff(7.701);
      d.getSmartNeighbours(from, to, I, J);
      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, I[0].size());
      ASSERT_EQ(1, I[0][0].size());
      EXPECT_FLOAT_EQ(1, I[0][0][0]);
      EXPECT_FLOAT_EQ(1, J[0][0][0]);

      // Use best neighbours (at 20m and 15m)
      d.setMinElevDiff(7.6999);
      d.getSmartNeighbours(from, to, I, J);
      ASSERT_EQ(1, I[0].size());
      ASSERT_EQ(2, I[0][0].size());
      EXPECT_FLOAT_EQ(2, I[0][0][0]);
      EXPECT_FLOAT_EQ(0, J[0][0][0]);
      EXPECT_FLOAT_EQ(0, I[0][0][1]);
      EXPECT_FLOAT_EQ(1, J[0][0][1]);

   }
   TEST_F(TestDownscalerSmart, fewerPointsThanSmart) {
      FileFake from(Options("nLat=5 nLon=3 nEns=1 nTime=1"));
      setLatLonElev(from, (const float[]) {4,5,6,7,8}, (const float[]){5,10,15}, (const float[]){70,50,20,80,60,70,50,40,30,20,10,40,50,30,60});
      FileFake to(Options("nLat=1 nLon=4 nEns=1 nTime=1"));
      setLatLonElev(to, (const float[]) {5.5}, (const float[]){2,4, 10,20}, (const float[]){120, 80, 600, 600});
      vec3Int I, J;

      DownscalerSmart d(Variable::Precip, Options());
      d.setSearchRadius(1);
      d.setNumSmart(20);
      d.getSmartNeighbours(from, to, I, J);

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

      DownscalerSmart d(Variable::Precip, Options());
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

      DownscalerSmart d(Variable::Precip, Options());
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

      DownscalerSmart d(Variable::Precip, Options());
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
