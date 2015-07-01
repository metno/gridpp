#include "../ParameterFile/ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   TEST(ParameterFileNetcdfTest, singleTime) {
      ParameterFileNetcdf file(Options("file=testing/files/10x10_param.nc"));
      ASSERT_EQ(2, file.getTimes().size());
      // Location loc(5,5,142.1456);
      Location loc(5,5,3);
      Parameters par = file.getParameters(0, loc);
      ASSERT_EQ(2, par.size());
      EXPECT_FLOAT_EQ(0.3, par[0]);
      EXPECT_FLOAT_EQ(2.3, par[1]);
   }
   TEST(ParameterFileNetcdfTest, description) {
      ParameterFileNetcdf::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
