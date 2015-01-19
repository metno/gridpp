#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile.h"
#include "../Calibrator/Zaga.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorZaga : public ::testing::Test {
      protected:
         TestCalibratorZaga() {
         }
         virtual ~TestCalibratorZaga() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
         void test(CalibratorZaga& cal, File& file, float i0, float i1, float i2, float e0, float e1, float e2) {
            FieldPtr field  = file.getField(Variable::Precip, 0);
            (*field)[0][0][0] = i0;
            (*field)[0][0][1] = i1;
            (*field)[0][0][2] = i2;
            cal.calibrate(file);
            const FieldPtr after  = file.getField(Variable::Precip, 0);

            EXPECT_EQ(e0, (*after)[0][0][0]);
            EXPECT_EQ(e1, (*after)[0][0][1]);
            EXPECT_EQ(e2, (*after)[0][0][2]);
         }
         void testP0(float ensMean, float ensFrac, float par4, float par5, float par6, float par7, float expected) {
            std::vector<float> parValues(8, 0);
            parValues[4] = par4;
            parValues[5] = par5;
            parValues[6] = par6;
            parValues[7] = par7;
            Parameters par(parValues);
            float p0 = CalibratorZaga::getP0(ensMean, ensFrac, par);
            EXPECT_FLOAT_EQ(expected, p0);
         }
   };

   TEST_F(TestCalibratorZaga, small) {
      // Set up file
      FileFake file(1, 1, 3, 1);
      FieldPtr field = file.getField(Variable::Precip, 0);
      (*field)[0][0][0] = 0;
      (*field)[0][0][1] = 0;
      (*field)[0][0][2] = 1;
      const Field before = *file.getField(Variable::Precip, 0);

      // Set up parameters
      ParameterFile parFile("Testing/files/parameters.txt");
      std::vector<float> parValues(8, 0);
      parValues[0] = 0;
      parValues[1] = 1;
      parValues[2] = 0;
      parValues[3] = 1;
      parValues[4] = 0;
      parValues[5] = 0;
      parValues[6] = 0;
      parValues[7] = 0;
      Parameters par(parValues);
      parFile.setParameters(par, 0);
      CalibratorZaga cal(parFile, Variable::Precip);
      test(cal, file, 12, 1, 4, 2, 3, 1);
   }

   TEST_F(TestCalibratorZaga, p0) {
      testP0(5, 0, 0, 0, 0, 0, 0.5);
      testP0(5, 0.5, 0, 0, 0, 0, 0.5);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
