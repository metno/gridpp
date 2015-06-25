#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Accumulate.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorAccumlate : public ::testing::Test {
      protected:
         TestCalibratorAccumlate() {
         }
         virtual ~TestCalibratorAccumlate() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
   };

   TEST_F(TestCalibratorAccumlate, 1x1) {
      Variable::Type var = Variable::T;
      Variable::Type varAcc = Variable::T;
      FileArome from("testing/files/1x1.nc");
      CalibratorAccumulate cal = CalibratorAccumulate(var, Options("outputVariable=T"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(0, (*from.getField(varAcc, 0))(0,0,0));
      EXPECT_FLOAT_EQ(20, (*from.getField(varAcc, 1))(0,0,0));
      EXPECT_FLOAT_EQ(35, (*from.getField(varAcc, 2))(0,0,0));
      EXPECT_FLOAT_EQ(56, (*from.getField(varAcc, 3))(0,0,0));
      EXPECT_FLOAT_EQ(70, (*from.getField(varAcc, 4))(0,0,0));
      EXPECT_FLOAT_EQ(100, (*from.getField(varAcc, 5))(0,0,0));
      EXPECT_FLOAT_EQ(121, (*from.getField(varAcc, 6))(0,0,0));
      EXPECT_FLOAT_EQ(140, (*from.getField(varAcc, 7))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(varAcc, 8))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(varAcc, 9))(0,0,0));
   }
   TEST_F(TestCalibratorAccumlate, 10x10) {
      Variable::Type var = Variable::Precip;
      Variable::Type varAcc = Variable::PrecipAcc;
      FileArome from("testing/files/10x10.nc");
      CalibratorAccumulate cal = CalibratorAccumulate(var, Options(""));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(0, (*from.getField(varAcc, 0))(5,2,0));
      EXPECT_FLOAT_EQ(0.539526, (*from.getField(varAcc, 1))(5,2,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(varAcc, 0))(5,9,0));
      EXPECT_FLOAT_EQ(6.929162, (*from.getField(varAcc, 1))(5,9,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(varAcc, 0))(0,9,0));
      EXPECT_FLOAT_EQ(5.442121, (*from.getField(varAcc, 1))(0,9,0));
   }
   TEST_F(TestCalibratorAccumlate, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      FileArome from("testing/files/1x1.nc");
      // File does not contain QNH so it should not be possible to accumulate
      CalibratorAccumulate cal = CalibratorAccumulate(Variable::QNH, Options(""));
      EXPECT_DEATH(cal.calibrate(from), ".*");
      CalibratorAccumulate cal2 = CalibratorAccumulate(Variable::QNH, Options("outputVariable=T"));
      EXPECT_DEATH(cal2.calibrate(from), ".*");
   }
   TEST_F(TestCalibratorAccumlate, description) {
      CalibratorAccumulate::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
