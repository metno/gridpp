#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Diagnose.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorDiagnose : public ::testing::Test {
      protected:
         TestCalibratorDiagnose() {
         }
         virtual ~TestCalibratorDiagnose() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
   };
   TEST_F(TestCalibratorDiagnose, WindSpeed) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorDiagnose cal(Variable::W ,Options(""));

      cal.calibrate(from, NULL);

      FieldPtr after = from.getField(Variable::W, 0);
      EXPECT_FLOAT_EQ(1.889807, (*after)(5,2,0)); // x: -1.886 y: -0.108
   }
   TEST_F(TestCalibratorDiagnose, U_V) {
      // Diagnose windspeed and direction from u and v. Then change the windspeed to
      // 13 and check that rediagnosis of U/V yields correct winds
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorDiagnose calW(Variable::W,Options(""));
      CalibratorDiagnose calWD(Variable::WD,Options(""));
      CalibratorDiagnose calU(Variable::U,Options(""));
      CalibratorDiagnose calV(Variable::V,Options(""));
      CalibratorDiagnose calXwind(Variable::Xwind,Options(""));
      CalibratorDiagnose calYwind(Variable::Ywind,Options(""));

      calW.calibrate(from, NULL);
      calWD.calibrate(from, NULL);

      // Diagnose W
      FieldPtr fieldW = from.getField(Variable::W, 0);
      EXPECT_FLOAT_EQ(1.8898070623831842, (*fieldW)(5,2,0));
      (*fieldW)(5,2,0) = 13;

      // Re-diagnose the new U and V
      calU.calibrate(from, NULL);
      calV.calibrate(from, NULL);
      calXwind.calibrate(from, NULL);
      calYwind.calibrate(from, NULL);

      FieldPtr fieldU = from.getField(Variable::U, 0);
      FieldPtr fieldV = from.getField(Variable::V, 0);
      FieldPtr fieldXwind = from.getField(Variable::Xwind, 0);
      FieldPtr fieldYwind = from.getField(Variable::Ywind, 0);
      EXPECT_FLOAT_EQ(-12.978647, (*fieldU)(5,2,0));
      EXPECT_FLOAT_EQ(-0.74479485, (*fieldV)(5,2,0));
      EXPECT_FLOAT_EQ(-12.978647, (*fieldXwind)(5,2,0));
      EXPECT_FLOAT_EQ(-0.74479485, (*fieldYwind)(5,2,0));
   }
   TEST_F(TestCalibratorDiagnose, dewpoint2RH) {
      // EXPECT_FLOAT_EQ(0.9, CalibratorDiagnose::dewpoint2RH(273.15, 273.15-1.27));
      EXPECT_FLOAT_EQ(1, CalibratorDiagnose::dewpoint2RH(273.15, 273.15));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnose::dewpoint2RH(Util::MV, 273.15-1.27));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnose::dewpoint2RH(273.15, Util::MV));
   }
   TEST_F(TestCalibratorDiagnose, RH2dewpoint) {
      // EXPECT_FLOAT_EQ(273.15-1.27, CalibratorDiagnose::RH2dewpoint(273.15, 0.9));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnose::RH2dewpoint(Util::MV, 0.9));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnose::RH2dewpoint(273.15, Util::MV));
   }
   TEST_F(TestCalibratorDiagnose, both) {
      // EXPECT_FLOAT_EQ(0.9, CalibratorDiagnose::dewpoint2RH(273.15, CalibratorDiagnose::RH2dewpoint(273.15, 0.9)));
      // EXPECT_FLOAT_EQ(1, CalibratorDiagnose::dewpoint2RH(273.15, CalibratorDiagnose::RH2dewpoint(273.15, 1)));
   }
   TEST_F(TestCalibratorDiagnose, DontKnowHow) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorDiagnose cal(Variable::T ,Options(""));
      EXPECT_DEATH(cal.calibrate(from, NULL), ".*");
   }
   TEST_F(TestCalibratorDiagnose, description) {
      CalibratorDiagnose::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
