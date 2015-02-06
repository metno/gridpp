#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile.h"
#include "../Parameters.h"
#include "../Calibrator/Smooth.h"
#include <gtest/gtest.h>
#include <cmath>

namespace {
   class TestCalibratorWindDirection : public ::testing::Test {
      protected:
         TestCalibratorWindDirection() {
         }
         virtual ~TestCalibratorWindDirection() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
         Parameters getParameters(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8, float a9) {
            std::vector<float> parValues(9, 0);
            parValues[0] = a1;
            parValues[1] = a2;
            parValues[2] = a3;
            parValues[3] = a4;
            parValues[4] = a5;
            parValues[5] = a6;
            parValues[6] = a7;
            parValues[7] = a8;
            parValues[8] = a9;
            return Parameters (parValues);
         }
         ParameterFile getParameterFile(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8, float a9) {
            ParameterFile parFile("testing/files/parametersWindDirection.txt");
            parFile.setParameters(getParameters(a1, a2, a3, a4, a5, a6, a7, a8, a9), 0);
            return parFile;
         }
         CalibratorWindDirection getCalibrator(ParameterFile* parFile) {
            return CalibratorWindDirection(parFile, Variable::T);
         }
   };
   /*
   TEST_F(TestCalibratorWindDirection, 10x10_noCorrection) {
      // Test that when the correction factor is 1, the values are unchanged
      FileArome from("testing/files/10x10.nc");

      ParameterFile parFile = getParameterFile(1,0,0,0,0,0,0,0,0);
      CalibratorWindDirection cal = getCalibrator(&parFile);

      Field before = *from.getField(Variable::T, 0);
      cal.calibrate(from);
      Field after   = *from.getField(Variable::T, 0);
      ASSERT_EQ(10, after.getNumLat());
      ASSERT_EQ(10, after.getNumLon());
      ASSERT_EQ(1,  after.getNumEns());

      for(int i = 0; i < from.getNumLat(); i++) {
         for(int j = 0; j < from.getNumLon(); j++) {
            EXPECT_FLOAT_EQ(before(i,j,0), after(i,j,0));
         }
      }
   }
   TEST_F(TestCalibratorWindDirection, 10x10) {
      FileArome from("testing/files/10x10.nc");

      ParameterFile parFile = getParameterFile(1,0.5,0.2,0,0,0,0,0,0);
      CalibratorWindDirection cal = getCalibrator(&parFile);

      Field before = *from.getField(Variable::T, 0);
      cal.calibrate(from);
      Field after   = *from.getField(Variable::T, 0);
      ASSERT_EQ(10, after.getNumLat());
      ASSERT_EQ(10, after.getNumLon());
      ASSERT_EQ(1,  after.getNumEns());

      for(int i = 0; i < from.getNumLat(); i++) {
         for(int j = 0; j < from.getNumLon(); j++) {
            EXPECT_FLOAT_EQ(before(i,j,0), after(i,j,0));
         }
      }
   }
   TEST_F(TestCalibratorWindDirection, 10x10_missingValues) {
      FileArome from("testing/files/10x10.nc");
      CalibratorSmooth cal = CalibratorSmooth(Variable::T);
      cal.setSmoothRadius(1);
      FieldPtr field = from.getField(Variable::T, 0);
      (*field)(4,1,0) = Util::MV;
      (*field)(4,2,0) = Util::MV;
      (*field)(4,3,0) = Util::MV;
      (*field)(5,1,0) = Util::MV;
      (*field)(5,2,0) = Util::MV;
      (*field)(5,3,0) = Util::MV;
      (*field)(6,1,0) = Util::MV;
      (*field)(6,2,0) = Util::MV;
      (*field)(6,3,0) = Util::MV;

      cal.calibrate(from);
      EXPECT_FLOAT_EQ(Util::MV, (*field)(5,2,0));
      EXPECT_FLOAT_EQ(304.6667, (*field)(5,1,0));
   }
   */
   TEST_F(TestCalibratorWindDirection, getFactor) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      ParameterFile parFile = getParameterFile(1,0,0,0,0,0,0,0,0);
      CalibratorWindDirection cal = getCalibrator(&parFile);
      Parameters par = parFile.getParameters(0);

      EXPECT_FLOAT_EQ(1, cal.getFactor(0, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(45, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(90, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(180, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(270, par));

      parFile = getParameterFile(1,0,1,0,0,0,0,0,0);
      cal = getCalibrator(&parFile);
      par = parFile.getParameters(0);

      EXPECT_FLOAT_EQ(2, cal.getFactor(0, par));
      EXPECT_FLOAT_EQ(1.7071067, cal.getFactor(45, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(90, par));
      EXPECT_NEAR(0, cal.getFactor(180, par), 1e-6);
      EXPECT_FLOAT_EQ(1, cal.getFactor(270, par));

      parFile = getParameterFile(1,1,0.5,0,0,0,0,0,0);
      cal = getCalibrator(&parFile);
      par = parFile.getParameters(0);

      EXPECT_FLOAT_EQ(1.5, cal.getFactor(0, par));
      EXPECT_FLOAT_EQ(2.06066, cal.getFactor(45, par));
      EXPECT_FLOAT_EQ(2, cal.getFactor(90, par));
      EXPECT_FLOAT_EQ(0.5, cal.getFactor(180, par));
      EXPECT_NEAR(0, cal.getFactor(270, par), 1e-6);
   }
   TEST_F(TestCalibratorWindDirection, getFactorInvalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      // Test that getFactor does not give negative values
      ParameterFile parFile = getParameterFile(1,2,3,0,0,0,0,0,0);
      CalibratorWindDirection cal = getCalibrator(&parFile);
      Parameters par = parFile.getParameters(0);
      for(int dir = 0; dir <= 360; dir = dir + 15) {
         EXPECT_GE(cal.getFactor(dir, par), 0);
      }

      // Test that getFactor returns invalid factor when parameters are invalid
      parFile = getParameterFile(1,2,3,0,0,0,Util::MV,0,0);
      cal = getCalibrator(&parFile);
      par = parFile.getParameters(0);
      for(int dir = 0; dir <= 360; dir = dir + 15) {
         EXPECT_FLOAT_EQ(Util::MV, cal.getFactor(dir, par));
      }

      // Invalid values
      // EXPECT_DEATH(cal.getFactor(-1, par), ".*");
   }
   TEST_F(TestCalibratorWindDirection, description) {
      CalibratorSmooth::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
