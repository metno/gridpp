#include "../ParameterFile/ParameterFile.h"
#include "../KalmanFilter.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   class KalmanFilterTest : public ::testing::Test {
      protected:
         KalmanFilterTest() {
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
   };

   TEST_F(KalmanFilterTest, default) {
      KalmanFilter kf(Variable::T, Options());
      ParameterFileText data(Options("file=testing/files/sampleObsFcst.txt"));
   }
   TEST_F(KalmanFilterTest, simple) {
      KalmanFilter kfSlow(Variable::T, Options("ratio=0.1"));
      KalmanFilter kfFast(Variable::T, Options("ratio=0.2"));
      float bias = 3;
      int   time = 0;
      KalmanParameters par0 = kfFast.initialize();
      KalmanParameters parFast = kfFast.update(bias, time, par0);
      KalmanParameters parSlow = kfSlow.update(bias, time, par0);
      ASSERT_EQ(par0.x.size(), parFast.x.size());
      ASSERT_EQ(par0.k.size(), parFast.k.size());
      ASSERT_EQ(par0.P.size(), parFast.P.size());
      ASSERT_EQ(par0.x.size(), parSlow.x.size());
      ASSERT_EQ(par0.k.size(), parSlow.k.size());
      ASSERT_EQ(par0.P.size(), parSlow.P.size());
      for(int i = 0; i < parFast.x.size(); i++) {
         // The estimated bias should be less than the measured bias
         EXPECT_LT(parFast.x[i], 3);
         EXPECT_GT(parFast.x[i], 0);
         EXPECT_LT(parSlow.x[i], 3);
         EXPECT_GT(parSlow.x[i], 0);
         // Kalman gain should be between 0 and 1
         EXPECT_LT(parFast.k[i], 1);
         EXPECT_GT(parFast.k[i], 0);
         EXPECT_LT(parSlow.k[i], 1);
         EXPECT_GT(parSlow.k[i], 0);
         // The fast KF should adapt quicker to the new bias
         EXPECT_LT(parSlow.x[i], parFast.x[i]);
         EXPECT_LT(parSlow.k[i], parFast.k[i]);
      }
   }
   /*
   TEST_F(KalmanFilterTest, calcBias) {
      KalmanFilter kf(Variable::T, 0.1);
      Parameters par = getParameters(0,1,2,3,4,5,6,7);
      kf.calcBias(par);
   }
   TEST_F(KalmanFilterTest, update) {
      KalmanFilter kf(Variable::T, 0.1);
      Parameters oldPar = getParameters(0,1,2,3,4,5,6,7);
      float obs = 3;
      float fcst = 2;
      Parameters newPar = kf.update(obs, fcst, oldPar);
      EXPECT_EQ(8, newPar.size());
   }
   */
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
