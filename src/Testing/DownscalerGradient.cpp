#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

namespace {
   class TestDownscalerGradient : public ::testing::Test {
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

   TEST_F(TestDownscalerGradient, 10x10) {
      DownscalerGradient d(Variable::T);
      d.setSearchRadius(1);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      FileFake to(1,4,1,from.getNumTime());
      setLatLonElev(to, (float[]) {5}, (float[]){2,2,12,20}, (float[]){120, 1500, 600, -100});
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT.getNumLat());
      ASSERT_EQ(4, toT.getNumLon());
      // T = T(nn) + gradient * (elev - elev(nn))
      EXPECT_FLOAT_EQ(301.31491, toT(0,0,0)); // 301 - 0.00797 * (120-160)
      EXPECT_FLOAT_EQ(290.34964, toT(0,1,0));
      EXPECT_FLOAT_EQ(301.29544, toT(0,2,0)); // 301 - 0.01068 * (600-346)
      EXPECT_FLOAT_EQ(308.77686, toT(0,3,0));

      // Fix the gradient
      d.setConstantGradient(0.01);
      d.downscale(from, to);
      const Field& toT2  = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT2.getNumLat());
      ASSERT_EQ(4, toT2.getNumLon());
      // T = T(nn) + gradient * (elev - elev(nn))
      EXPECT_FLOAT_EQ(300.60367, toT2(0,0,0));
      EXPECT_FLOAT_EQ(314.40369, toT2(0,1,0));
      EXPECT_FLOAT_EQ(306.53052, toT2(0,2,0));
      EXPECT_FLOAT_EQ(299.53052, toT2(0,3,0));
   }
   // Check that logTransform works
   TEST_F(TestDownscalerGradient, 10x10log) {
      DownscalerGradient d(Variable::T);
      d.setSearchRadius(1);
      d.setLogTransform(true);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      FileFake to(1,4,1,from.getNumTime());
      setLatLonElev(to, (float[]) {5}, (float[]){2,2,12,20}, (float[]){120, 1500, 600, -100});
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT.getNumLat());
      ASSERT_EQ(4, toT.getNumLon());
      // T = T(nn) * exp(gradient * (elev - elev(nn)))
      EXPECT_FLOAT_EQ(301.30985, toT(0,0,0)); // 301 * exp(-2.59620e-5 * (120-159.63))

      // Fix the gradient
      d.setConstantGradient(-0.01);
      d.setLogTransform(true);
      d.downscale(from, to);
      const Field& toT2  = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT2.getNumLat());
      ASSERT_EQ(4, toT2.getNumLon());
      // T = T(nn) * exp(gradient * (elev - elev(nn)))
      EXPECT_FLOAT_EQ(447.3916, toT2(0,0,0));
   }
   TEST_F(TestDownscalerGradient, 10x10minElevDiff) {
      DownscalerGradient d(Variable::T);
      d.setSearchRadius(1);
      d.setMinElevDiff(1500);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      FileFake to(1,4,1,from.getNumTime());
      setLatLonElev(to, (float[]) {5}, (float[]){2,2,12,20}, (float[]){120, 1500, 600, -100});
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT.getNumLat());
      ASSERT_EQ(4, toT.getNumLon());
      // T = T(nn) + gradient * (elev - elev(nn))
      EXPECT_FLOAT_EQ(301, toT(0,0,0));
      EXPECT_FLOAT_EQ(301, toT(0,1,0));
      EXPECT_FLOAT_EQ(304, toT(0,2,0));
      EXPECT_FLOAT_EQ(304, toT(0,3,0));

      // Fix the gradient
      // This should not be affected by the minElevDiff
      d.setConstantGradient(0.01);
      d.downscale(from, to);
      const Field& toT2  = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT2.getNumLat());
      ASSERT_EQ(4, toT2.getNumLon());
      // T = T(nn) + gradient * (elev - elev(nn))
      EXPECT_FLOAT_EQ(300.60367, toT2(0,0,0));
      EXPECT_FLOAT_EQ(314.40369, toT2(0,1,0));
      EXPECT_FLOAT_EQ(306.53052, toT2(0,2,0));
      EXPECT_FLOAT_EQ(299.53052, toT2(0,3,0));
   }
   TEST_F(TestDownscalerGradient, 10x10negativeTemperatures) {
      DownscalerGradient d(Variable::T);
      d.setSearchRadius(1);
      d.setMinElevDiff(0);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      FileFake to(1,3,1,from.getNumTime());
      setLatLonElev(to, (float[]) {5}, (float[]){2,2,2}, (float[]){100000, 10000, 0});
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT.getNumLat());
      ASSERT_EQ(3, toT.getNumLon());
      // Gradient = -0.00797
      EXPECT_FLOAT_EQ(301,   toT(0,0,0)); // nearest neighbour
      EXPECT_FLOAT_EQ(222.80988, toT(0,1,0));
      EXPECT_FLOAT_EQ(302.2684, toT(0,2,0));

      // Fix the gradient
      d.setConstantGradient(-0.1);
      d.downscale(from, to);
      const Field& toT2  = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT2.getNumLat());
      ASSERT_EQ(3, toT2.getNumLon());
      EXPECT_FLOAT_EQ(301, toT2(0,0,0)); // nearest neighbour
      EXPECT_FLOAT_EQ(301, toT2(0,1,0)); // nearest neighbour
      EXPECT_FLOAT_EQ(316.96323, toT2(0,2,0));
   }
   TEST_F(TestDownscalerGradient, missingValues) {

   }
   TEST_F(TestDownscalerGradient, constantGradient) {
      DownscalerGradient d(Variable::T);
      d.setSearchRadius(1);
      FileArome from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(Variable::T, 0);
      FileFake to(1,4,1,from.getNumTime());
      setLatLonElev(to, (float[]) {5}, (float[]){2,2,12,20}, (float[]){120, 1500, 600, -100});
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(Variable::T, 0);
      ASSERT_EQ(1, toT.getNumLat());
      ASSERT_EQ(4, toT.getNumLon());
      // T = T(nn) + gradient * (elev - elev(nn))
      EXPECT_FLOAT_EQ(301.31491, toT(0,0,0));
      EXPECT_FLOAT_EQ(290.34964,  toT(0,1,0));
      EXPECT_FLOAT_EQ(301.29544,  toT(0,2,0));
      EXPECT_FLOAT_EQ(308.77686,  toT(0,3,0));
   }
   TEST_F(TestDownscalerGradient, setGetSearchRadius) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      DownscalerGradient d(Variable::Precip);
      // Check that default is valid
      EXPECT_GE(d.getSearchRadius(), 0);
      d.setSearchRadius(5);
      EXPECT_FLOAT_EQ(5, d.getSearchRadius());

      // Invalid values
      EXPECT_DEATH(d.setSearchRadius(-1), ".*");
      EXPECT_DEATH(d.setSearchRadius(0), ".*");
      EXPECT_DEATH(d.setSearchRadius(Util::MV), ".*");
   }
   TEST_F(TestDownscalerGradient, setGetConstantGradient) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      DownscalerGradient d(Variable::Precip);
      d.setConstantGradient(5);
      EXPECT_FLOAT_EQ(5, d.getConstantGradient());
      d.setConstantGradient(1);
      EXPECT_FLOAT_EQ(1, d.getConstantGradient());

      // Invalid values
      EXPECT_DEATH(d.setConstantGradient(Util::MV), ".*");
   }
   TEST_F(TestDownscalerGradient, setGetMinElevDiff) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      DownscalerGradient d(Variable::Precip);
      // Check that default is valid
      EXPECT_GE(d.getMinElevDiff(), 0);
      d.setMinElevDiff(213.1);
      EXPECT_FLOAT_EQ(213.1, d.getMinElevDiff());
   }
   TEST_F(TestDownscalerGradient, description) {
      DownscalerGradient::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
