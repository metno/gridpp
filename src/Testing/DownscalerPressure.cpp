#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Pressure.h"
#include <gtest/gtest.h>
#include "gridpp.h"

namespace {
   class TestDownscalerPressure : public ::testing::Test {
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
            int nLat = iFile.getNumY();
            int nLon = iFile.getNumX();
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
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(TestDownscalerPressure, 10x10) {
      DownscalerPressure d(mVariable, mVariable, Options());
      FileNetcdf from("testing/files/10x10.nc");
      std::stringstream ss;
      ss << "nLat=1 nLon=4 nEns=1 nTime=" << from.getNumTime();
      FileFake to(Options(ss.str()));
      setLatLonElev(to, (const float[]) {5}, (const float[]){2,2,12,20}, (const float[]){120, 1500, 600, -100});
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(mVariable, 0);
      ASSERT_EQ(1, toT.getNumY());
      ASSERT_EQ(4, toT.getNumX());
      // T = T(nn) * exp(-1.21e-4 * (elev - elev(nn)))
      EXPECT_FLOAT_EQ(302.41766, toT(0,0,0)); // 301 160m->120m
      EXPECT_FLOAT_EQ(256.77454, toT(0,1,0)); // 301 160m->1500m
      EXPECT_FLOAT_EQ(295.01501, toT(0,2,0)); // 304 347m->600m
      EXPECT_FLOAT_EQ(320.54321, toT(0,3,0)); // 304 347m->-100m
   }
   TEST_F(TestDownscalerPressure, calcPressure) {
      EXPECT_FLOAT_EQ(101325, gridpp::pressure(0, 0, 101325));
      EXPECT_FLOAT_EQ(89996.852, gridpp::pressure(0, 1000, 101325));
      EXPECT_FLOAT_EQ(67554.031, gridpp::pressure(300, 600, 70000));

      EXPECT_FLOAT_EQ(Util::MV, gridpp::pressure(gridpp::MV, 100, 101325));
      EXPECT_FLOAT_EQ(Util::MV, gridpp::pressure(0, 100, gridpp::MV));
      EXPECT_FLOAT_EQ(Util::MV, gridpp::pressure(0, gridpp::MV, 101325));
   }
   TEST_F(TestDownscalerPressure, description) {
      DownscalerPressure::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
