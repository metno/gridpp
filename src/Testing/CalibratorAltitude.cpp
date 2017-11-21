#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Altitude.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorAltitude : public ::testing::Test {
      protected:
         TestCalibratorAltitude() {
         }
         virtual ~TestCalibratorAltitude() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
   };
   TEST_F(TestCalibratorAltitude, arome) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileNetcdf par(Options("file=testing/files/10x10_param_zero_altitude.nc"));
      CalibratorAltitude cal = CalibratorAltitude(Options());

      cal.calibrate(from, &par);

      // Elevations should all be 0
      vec2 elevs = from.getElevs();
      for(int i = 0; i < from.getNumLat(); i++) {
         for(int j = 0; j < from.getNumLon(); j++) {
            EXPECT_FLOAT_EQ(0, elevs[i][j]);
         }
      }
      // Shouldn't have changed anything else
      FieldPtr after = from.getField(Variable::T, 0);
      EXPECT_FLOAT_EQ(301, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(304, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(320, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorAltitude, ec) {
      FileEc from("testing/files/10x10_ec.nc");
      ParameterFileNetcdf par(Options("file=testing/files/10x10_param_zero_altitude.nc"));
      CalibratorAltitude cal = CalibratorAltitude(Options());

      cal.calibrate(from, &par);

      // Elevations should all be 0
      vec2 elevs = from.getElevs();
      for(int i = 0; i < from.getNumLat(); i++) {
         for(int j = 0; j < from.getNumLon(); j++) {
            EXPECT_FLOAT_EQ(0, elevs[i][j]);
         }
      }
      // Shouldn't have changed anything else
      FieldPtr after = from.getField(Variable::T, 0);
      EXPECT_FLOAT_EQ(301, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(304, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(320, (*after)(0,9,0));
   }
   // Should not calibrate when a parameter file with no locaions is used
   TEST_F(TestCalibratorAltitude, locationIndependent) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileText par(Options("file=testing/files/parametersSingleTime.txt"));
      CalibratorAltitude cal = CalibratorAltitude(Options());
      EXPECT_DEATH(cal.calibrate(from, &par), ".*");
   }
   TEST_F(TestCalibratorAltitude, description) {
      CalibratorAltitude::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
