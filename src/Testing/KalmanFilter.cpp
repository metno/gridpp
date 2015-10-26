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
      ParameterFileText data(Options("file=testing/files/sampleObsFcst.txt spatial=1"));
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
