#define BOOST_TEST_MODULE Base
#include <gtest/gtest.h>
#include <math.h>
#include "Calibration.h"
#include "ParameterFile.h"

namespace {
   class BaseTest : public ::testing::Test {
      protected:
         BaseTest() {
            // You can do set-up work for each test here.
         }
         virtual ~BaseTest() {
            // You can do clean-up work that doesn't throw exceptions here.
         }
         virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
         }

         virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
         }

   };

   TEST(Calibration, logit) {
      const float p[] = {0.2,0.9};
      const float exp[] = {-1.386294, 2.197225};
      for(int i = 0; i < sizeof(p)/sizeof(float); i++) {
         float ans = Calibration::logit(p[i]);
         EXPECT_FLOAT_EQ(exp[i], ans);
      }
   }
   TEST(Calibration, invlogit) {
      const float x[] = {-2, 3};
      const float exp[] = {0.1192029, 0.9525741};
      for(int i = 0; i < sizeof(x)/sizeof(float); i++) {
         float ans = Calibration::invLogit(x[i]);
         EXPECT_FLOAT_EQ(exp[i], ans);
      }
   }
   TEST(ParameterFile, singleTime) {
      ParameterFile file("testFiles/parameterSingleTime.txt");
      const float x[] = {-2, 3};
      const float exp[] = {0.1192029, 0.9525741};
      for(int i = 0; i < sizeof(x)/sizeof(float); i++) {
         float ans = Calibration::invLogit(x[i]);
         EXPECT_FLOAT_EQ(exp[i], ans);
      }
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
