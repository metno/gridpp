#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/DiagnoseWindSpeed.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorDiagnoseWindSpeed : public ::testing::Test {
      protected:
         virtual void SetUp() {
            mX = Variable("x_wind_10m");
            mY = Variable("y_wind_10m");
            mW = Variable("windspeed_10m");
            cal = new CalibratorDiagnoseWindSpeed(mW ,Options("x=x_wind_10m y=y_wind_10m"));
         }
         virtual void TearDown() {
            delete cal;
         }
         Variable mX;
         Variable mY;
         Variable mW;
         CalibratorDiagnoseWindSpeed* cal;
   };
   TEST_F(TestCalibratorDiagnoseWindSpeed, U_V) {
      FileNetcdf from("testing/files/10x10.nc");
      FieldPtr fieldX = from.getField(mX, 0);
      FieldPtr fieldY = from.getField(mY, 0);
      // Set some values to missing
      (*fieldX)(1, 2, 0) = Util::MV;
      (*fieldY)(2, 2, 0) = Util::MV;
      (*fieldY)(3, 2, 0) = Util::MV;

      cal->calibrate(from, NULL);
      FieldPtr fieldW = from.getField(mW, 0);
      EXPECT_FLOAT_EQ(1.8898070623831842, (*fieldW)(5,2,0));  // x: -1.886 y: -0.108
      EXPECT_FLOAT_EQ(Util::MV, (*fieldW)(1,2,0));
      EXPECT_FLOAT_EQ(Util::MV, (*fieldW)(2,2,0));
      EXPECT_FLOAT_EQ(Util::MV, (*fieldW)(3,2,0));
   }
   TEST_F(TestCalibratorDiagnoseWindSpeed, description) {
      CalibratorDiagnoseWindSpeed::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
