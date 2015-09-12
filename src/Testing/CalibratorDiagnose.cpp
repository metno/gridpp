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
      FileArome from("testing/files/10x10.nc");
      CalibratorDiagnose cal(Variable::W ,Options(""));

      cal.calibrate(from, NULL);

      FieldPtr after = from.getField(Variable::W, 0);
      EXPECT_FLOAT_EQ(1.889807, (*after)(5,2,0)); // x: -1.886 y: -0.108
   }
   TEST_F(TestCalibratorDiagnose, CantDiagnose) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      FileArome from("testing/files/10x10.nc");
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
