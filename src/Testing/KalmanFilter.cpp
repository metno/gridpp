#include "../ParameterFile/ParameterFile.h"
#include "../KalmanFilter.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   TEST(KalmanFilter, default) {
      KalmanFilter kf(Variable::T, 0.1);
      ParameterFileText  data(Options("file=testing/files/sampleObsFcst.txt spatial=1"));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
