#include "../ParameterFile/ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   TEST(ParameterFileMetnoKalman, default) {
      ParameterFileMetnoKalman file("testing/files/kalmanOutput.txt");
      std::vector<Location> locations = file.getLocations();
      ASSERT_EQ(805, locations.size());
      Location locFirst = Location(70.9394,-8.6690,10);

      Parameters par = file.getParameters(0, locFirst);
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(0.48, par[0]);
      par = file.getParameters(3, locFirst);
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(0.51, par[0]);
      par = file.getParameters(66, locFirst);
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(0.49, par[0]);

      // Interpolated between (0, 0.48) and (3, 0.51)
      par = file.getParameters(1, locFirst);
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(0.49, par[0]);
      par = file.getParameters(2, locFirst);
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(0.50, par[0]);

      Location locLast(52.2000,27.8700,136);
      par = file.getParameters(0, locLast);
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(-1.59, par[0]);
      par = file.getParameters(66, locLast);
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(-1.05, par[0]);
   }
   TEST(ParameterFileMetnoKalman, emptyFile) {
      ParameterFileMetnoKalman file("testing/files/89h9382he9823he92.txt");
      std::vector<Location> locations = file.getLocations();
      EXPECT_EQ(0, locations.size());
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
