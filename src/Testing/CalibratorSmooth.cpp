#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile.h"
#include "../Calibrator/Smooth.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorSmooth : public ::testing::Test {
      protected:
         TestCalibratorSmooth() {
         }
         virtual ~TestCalibratorSmooth() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
   };
   TEST_F(TestCalibratorSmooth, 10x10_doubleSmooth) {
      FileArome from("testing/files/10x10.nc");
      CalibratorSmooth cal = CalibratorSmooth(Variable::T);
      cal.setSmoothRadius(1);

      cal.calibrate(from);
      FieldPtr after = from.getField(Variable::T, 0);
      ASSERT_EQ(10, (*after).size());
      ASSERT_EQ(10, (*after)[0].size());
      ASSERT_EQ(1,  (*after)[0][0].size());

      EXPECT_FLOAT_EQ(304.6667, (*after)[5][2][0]);
      EXPECT_FLOAT_EQ(306.1667, (*after)[5][9][0]);
      EXPECT_FLOAT_EQ(308.25,   (*after)[0][9][0]);

      cal.setSmoothRadius(2);
      cal.calibrate(from);
      EXPECT_FLOAT_EQ(304.73114, (*after)[5][2][0]);
      EXPECT_FLOAT_EQ(305.35556, (*after)[5][9][0]);
   }
   TEST_F(TestCalibratorSmooth, 10x10_missingValues) {
      FileArome from("testing/files/10x10.nc");
      CalibratorSmooth cal = CalibratorSmooth(Variable::T);
      cal.setSmoothRadius(1);
      FieldPtr field = from.getField(Variable::T, 0);
      (*field)[4][1][0] = Util::MV;
      (*field)[4][2][0] = Util::MV;
      (*field)[4][3][0] = Util::MV;
      (*field)[5][1][0] = Util::MV;
      (*field)[5][2][0] = Util::MV;
      (*field)[5][3][0] = Util::MV;
      (*field)[6][1][0] = Util::MV;
      (*field)[6][2][0] = Util::MV;
      (*field)[6][3][0] = Util::MV;

      cal.calibrate(from);
      EXPECT_FLOAT_EQ(Util::MV, (*field)[5][2][0]);
      EXPECT_FLOAT_EQ(304.6667, (*field)[5][1][0]);
   }
   TEST_F(TestCalibratorSmooth, setGetSmoothingRadius) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      CalibratorSmooth cal = CalibratorSmooth(Variable::T);
      // Check that default is valid
      EXPECT_GE(cal.getSmoothRadius(), 1);

      cal.setSmoothRadius(12);
      EXPECT_FLOAT_EQ(12, cal.getSmoothRadius());
      cal.setSmoothRadius(4);
      EXPECT_FLOAT_EQ(4, cal.getSmoothRadius());

      // Invalid values
      EXPECT_DEATH(cal.setSmoothRadius(-1), ".*");
      EXPECT_DEATH(cal.setSmoothRadius(0), ".*");
      EXPECT_DEATH(cal.setSmoothRadius(Util::MV), ".*");
   }
   TEST_F(TestCalibratorSmooth, description) {
      CalibratorSmooth::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
