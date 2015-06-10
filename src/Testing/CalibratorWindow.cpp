#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile.h"
#include "../Calibrator/Window.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorWindow : public ::testing::Test {
      protected:
         TestCalibratorWindow() {
         }
         virtual ~TestCalibratorWindow() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
   };
   TEST_F(TestCalibratorWindow, radius0) {
      FileArome from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(Variable::T ,Options("radius=0 operator=mean"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(23, (*from.getField(Variable::T, 0))(0,0,0));
      EXPECT_FLOAT_EQ(20, (*from.getField(Variable::T, 1))(0,0,0));
      EXPECT_FLOAT_EQ(15, (*from.getField(Variable::T, 2))(0,0,0));
      EXPECT_FLOAT_EQ(21, (*from.getField(Variable::T, 3))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(Variable::T, 4))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(Variable::T, 5))(0,0,0));
      EXPECT_FLOAT_EQ(21, (*from.getField(Variable::T, 6))(0,0,0));
      EXPECT_FLOAT_EQ(19, (*from.getField(Variable::T, 7))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(Variable::T, 8))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(Variable::T, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, radius2) {
      FileArome from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(Variable::T ,Options("radius=2 operator=mean"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(19.333333, (*from.getField(Variable::T, 0))(0,0,0));
      EXPECT_FLOAT_EQ(19.75,     (*from.getField(Variable::T, 1))(0,0,0));
      EXPECT_FLOAT_EQ(18.6,      (*from.getField(Variable::T, 2))(0,0,0));
      EXPECT_FLOAT_EQ(20,        (*from.getField(Variable::T, 3))(0,0,0));
      EXPECT_FLOAT_EQ(20.2,      (*from.getField(Variable::T, 4))(0,0,0));
      EXPECT_FLOAT_EQ(21,        (*from.getField(Variable::T, 5))(0,0,0));
      EXPECT_FLOAT_EQ(21,        (*from.getField(Variable::T, 6))(0,0,0));
      EXPECT_FLOAT_EQ(23.25,     (*from.getField(Variable::T, 7))(0,0,0));
      EXPECT_FLOAT_EQ(21,        (*from.getField(Variable::T, 8))(0,0,0));
      EXPECT_FLOAT_EQ(21,        (*from.getField(Variable::T, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, min) {
      FileArome from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(Variable::T ,Options("radius=2 operator=min"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(15, (*from.getField(Variable::T, 0))(0,0,0));
      EXPECT_FLOAT_EQ(15, (*from.getField(Variable::T, 1))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(Variable::T, 2))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(Variable::T, 3))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(Variable::T, 4))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(Variable::T, 5))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(Variable::T, 6))(0,0,0));
      EXPECT_FLOAT_EQ(19, (*from.getField(Variable::T, 7))(0,0,0));
      EXPECT_FLOAT_EQ(19, (*from.getField(Variable::T, 8))(0,0,0));
      EXPECT_FLOAT_EQ(19, (*from.getField(Variable::T, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, max) {
      FileArome from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(Variable::T ,Options("radius=2 operator=max"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(23, (*from.getField(Variable::T, 0))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(Variable::T, 1))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(Variable::T, 2))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(Variable::T, 3))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(Variable::T, 4))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(Variable::T, 5))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(Variable::T, 6))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(Variable::T, 7))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(Variable::T, 8))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(Variable::T, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, std) {
      FileArome from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(Variable::T ,Options("radius=2 operator=std"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(3.299832,  (*from.getField(Variable::T, 0))(0,0,0));
      EXPECT_FLOAT_EQ(2.947457,  (*from.getField(Variable::T, 1))(0,0,0));
      EXPECT_FLOAT_EQ(3.4985712, (*from.getField(Variable::T, 2))(0,0,0));
      EXPECT_FLOAT_EQ(4.145781,  (*from.getField(Variable::T, 7))(0,0,0));
      EXPECT_FLOAT_EQ(2,         (*from.getField(Variable::T, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, quantile) {
      FileArome from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(Variable::T ,Options("radius=2 operator=quantile quantile=0.5"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(20,   (*from.getField(Variable::T, 0))(0,0,0));
      EXPECT_FLOAT_EQ(20.5, (*from.getField(Variable::T, 1))(0,0,0));
      EXPECT_FLOAT_EQ(20,   (*from.getField(Variable::T, 2))(0,0,0));
      EXPECT_FLOAT_EQ(22,   (*from.getField(Variable::T, 7))(0,0,0));
      EXPECT_FLOAT_EQ(21,   (*from.getField(Variable::T, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, median) {
      FileArome from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(Variable::T ,Options("radius=2 operator=median"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(20,   (*from.getField(Variable::T, 0))(0,0,0));
      EXPECT_FLOAT_EQ(20.5, (*from.getField(Variable::T, 1))(0,0,0));
      EXPECT_FLOAT_EQ(20,   (*from.getField(Variable::T, 2))(0,0,0));
      EXPECT_FLOAT_EQ(22,   (*from.getField(Variable::T, 7))(0,0,0));
      EXPECT_FLOAT_EQ(21,   (*from.getField(Variable::T, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, radius100) {
      // A large radius forces all values to be the same
      FileArome from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(Variable::T ,Options("radius=100 operator=mean"));
      cal.calibrate(from);

      for(int t = 0; t < 10; t++) {
         EXPECT_FLOAT_EQ(20.666666, (*from.getField(Variable::T, t))(0,0,0));
      }
   }
   TEST_F(TestCalibratorWindow, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      EXPECT_DEATH(CalibratorWindow(Variable::T, Options("operator=invalidOperator")), ".*");

      // Negative radius
      EXPECT_DEATH(CalibratorWindow(Variable::T, Options("radius=-1")), ".*");
      EXPECT_DEATH(CalibratorWindow(Variable::T, Options("radius=-1 operator=quantile")), ".*");
      EXPECT_DEATH(CalibratorWindow(Variable::T, Options("radius=-1 operator=quantile quantile=0.5")), ".*");

      // Missing quantile
      EXPECT_DEATH(CalibratorWindow(Variable::T, Options("operator=quantile")), ".*");
      // Invalid quantile
      EXPECT_DEATH(CalibratorWindow(Variable::T, Options("operator=quantile quantile=-0.1")), ".*");
      EXPECT_DEATH(CalibratorWindow(Variable::T, Options("operator=quantile quantile=1.1")), ".*");
      EXPECT_DEATH(CalibratorWindow(Variable::T, Options("operator=quantile quantile=-999")), ".*");
   }
   TEST_F(TestCalibratorWindow, description) {
      CalibratorWindow::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
