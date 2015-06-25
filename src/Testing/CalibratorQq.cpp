#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Qq.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorQq : public ::testing::Test {
      protected:
         TestCalibratorQq() {
         }
         virtual ~TestCalibratorQq() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
         Parameters getParameters(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8) {
            std::vector<float> parValues(8, 0);
            parValues[0] = a1;
            parValues[1] = a2;
            parValues[2] = a3;
            parValues[3] = a4;
            parValues[4] = a5;
            parValues[5] = a6;
            parValues[6] = a7;
            parValues[7] = a8;
            return Parameters (parValues);
         }
         ParameterFile getParameterFile(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8) {
            ParameterFile parFile("testing/files/parameters.txt");
            parFile.setParameters(getParameters(a1, a2, a3, a4, a5, a6, a7, a8), 0);
            return parFile;
         }
   };
   TEST_F(TestCalibratorQq, 1to1) {
      FileArome from("testing/files/10x10.nc");
      ParameterFile parFile = getParameterFile(250,250,260,300,290,303,300,315);
      CalibratorQq cal(&parFile, Variable::T ,Options("extrapolation=1to1"));

      cal.calibrate(from);

      FieldPtr after = from.getField(Variable::T, 0);
      EXPECT_FLOAT_EQ(270, (*after)(5,2,0));      // 301 (row,col)
      EXPECT_FLOAT_EQ(290.8333, (*after)(5,9,0)); // 304
      EXPECT_FLOAT_EQ(305, (*after)(0,9,0));      // 320 // Outside
   }
   TEST_F(TestCalibratorQq, meanSlope) {
      FileArome from("testing/files/10x10.nc");
      ParameterFile parFile = getParameterFile(250,250,260,300,290,303,300,315);
      CalibratorQq cal(&parFile, Variable::T ,Options("extrapolation=meanSlope"));

      cal.calibrate(from);

      FieldPtr after = from.getField(Variable::T, 0);
      EXPECT_FLOAT_EQ(303.8462, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorQq, nearestSlope) {
      FileArome from("testing/files/10x10.nc");
      ParameterFile parFile = getParameterFile(250,250,260,300,290,303,300,315);
      CalibratorQq cal(&parFile, Variable::T ,Options("extrapolation=nearestSlope"));

      cal.calibrate(from);

      FieldPtr after = from.getField(Variable::T, 0);
      EXPECT_FLOAT_EQ(304.16666667, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorQq, twoEqual) {
      FileArome from("testing/files/10x10.nc");
      ParameterFile parFile = getParameterFile(250,250,270,301,280,301,320,320);
      CalibratorQq cal(&parFile, Variable::T ,Options("extrapolation=1to1"));

      cal.calibrate(from);

      FieldPtr after = from.getField(Variable::T, 0);
      EXPECT_FLOAT_EQ(275, (*after)(5,2,0));
   }
   TEST_F(TestCalibratorQq, twoEqualObs) {
      FileArome from("testing/files/10x10.nc");
      ParameterFile parFile = getParameterFile(250,250,280,300,290,302,320,320);
      CalibratorQq cal(&parFile, Variable::T ,Options("extrapolation=1to1"));

      cal.calibrate(from);

      FieldPtr after = from.getField(Variable::T, 0);
      EXPECT_FLOAT_EQ(285, (*after)(5,2,0)); // raw = 301
   }
   TEST_F(TestCalibratorQq, description) {
      CalibratorQc::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
