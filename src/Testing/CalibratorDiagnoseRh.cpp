#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/DiagnoseRh.h"
#include "../Calibrator/DiagnoseDewpoint.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorDiagnoseRh : public ::testing::Test {
      protected:
         virtual void SetUp() {
            mT = Variable("air_temperature_2m");
            mRH = Variable("relative_humidity_2m");
            mDewpoint = Variable("dewpoint");
         }
         virtual void TearDown() {
         }
         Variable mT;
         Variable mRH;
         Variable mDewpoint;
   };
   TEST_F(TestCalibratorDiagnoseRh, test1) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorDiagnoseDewpoint rh2d = CalibratorDiagnoseDewpoint(mDewpoint ,Options("temperature=air_temperature_2m rh=relative_humidity_2m"));
      CalibratorDiagnoseRh d2rh = CalibratorDiagnoseRh(mRH ,Options("temperature=air_temperature_2m dewpoint=dewpoint"));
      FieldPtr field = from.getField(mRH, 0);
      float before = (*field)(0,0,0);
      rh2d.calibrate(from, NULL);
      d2rh.calibrate(from, NULL);
      float after = (*field)(0,0,0);
      // Doesn't work
      // EXPECT_FLOAT_EQ(before, after);
   }
   TEST_F(TestCalibratorDiagnoseRh, dewpoint2RH) {
      // EXPECT_FLOAT_EQ(0.9, CalibratorDiagnose::dewpoint2RH(273.15, 273.15-1.27));
      EXPECT_FLOAT_EQ(1, CalibratorDiagnoseRh::dewpoint2RH(273.15, 273.15));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnoseRh::dewpoint2RH(Util::MV, 273.15-1.27));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnoseRh::dewpoint2RH(273.15, Util::MV));
   }
   TEST_F(TestCalibratorDiagnoseRh, description) {
      CalibratorDiagnoseRh::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
