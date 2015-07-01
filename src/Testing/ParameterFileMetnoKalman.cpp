#include "../ParameterFile/ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   TEST(ParameterFileMetnoKalman, default) {
      EXPECT_TRUE(ParameterFileMetnoKalman::isValid("testing/files/kalmanOutput.txt"));
      ParameterFileMetnoKalman file(Options("file=testing/files/kalmanOutput.txt"));
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

      EXPECT_EQ(67, file.getTimes().size());

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

      // Missing last obs value before biases (should not affect biases)
      par = file.getParameters(0, Location(71.0937, 23.9817, 0));
      EXPECT_FLOAT_EQ(0.15, par[0]);
   }
   TEST(ParameterFileMetnoKalman, missingValue) {
      EXPECT_TRUE(ParameterFileMetnoKalman::isValid("testing/files/kalmanOutputWithMissing.txt"));
      ParameterFileMetnoKalman file(Options("file=testing/files/kalmanOutputWithMissing.txt"));
      std::vector<Location> locations = file.getLocations();
      ASSERT_EQ(21, locations.size());
      Location loc = Location(69.0577,18.5437,0);

      Parameters par = file.getParameters(0, loc);
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(0.37, par[0]);

      // Missing value at 3 index, corresponding to time 6
      // Nearby interpolated values should also be missing
      for(int t = 4; t <= 8; t++) {
         par = file.getParameters(t, loc);
         ASSERT_EQ(1, par.size());
         EXPECT_FLOAT_EQ(Util::MV, par[0]);
      }
   }
   TEST(ParameterFileMetnoKalman, invalidFile) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // Wrongly formatted file
      EXPECT_FALSE(ParameterFileMetnoKalman::isValid("testing/files/parameters.txt"));
      EXPECT_DEATH(ParameterFileMetnoKalman(Options("testing/files/parameters.txt")), ".*");
      // Empty file
      EXPECT_FALSE(ParameterFileMetnoKalman::isValid("testing/files/nonexistaiowenwenrewoi.txt"));
      EXPECT_DEATH(ParameterFileMetnoKalman(Options("testing/files/nonexistaiowenwenrewoi.txt")), ".*");
      // One missing value on a row
      EXPECT_FALSE(ParameterFileMetnoKalman::isValid("testing/files/kalmanInvalid1.txt"));
      EXPECT_DEATH(ParameterFileMetnoKalman(Options("testing/files/kalmanInvalid1.txt")), ".*");
      // Wrong number of columns (header says 20, in reallity its 23)
      EXPECT_FALSE(ParameterFileMetnoKalman::isValid("testing/files/kalmanInvalid2.txt"));
      EXPECT_DEATH(ParameterFileMetnoKalman(Options("testing/files/kalmanInvalid2.txt")), ".*");
      // Missing number of times in header
      EXPECT_FALSE(ParameterFileMetnoKalman::isValid("testing/files/kalmanInvalid3.txt"));
      EXPECT_DEATH(ParameterFileMetnoKalman(Options("testing/files/kalmanInvalid3.txt")), ".*");
      // Text in station ID
      EXPECT_FALSE(ParameterFileMetnoKalman::isValid("testing/files/kalmanInvalid4.txt"));
      EXPECT_DEATH(ParameterFileMetnoKalman(Options("testing/files/kalmanInvalid4.txt")), ".*");
   }
   TEST(ParameterFileMetnoKalman, emptyFile) {
      ParameterFileMetnoKalman file(Options("file=testing/files/kalmanEmpty.txt"));
      std::vector<Location> locations = file.getLocations();
      EXPECT_EQ(0, locations.size());
   }
   TEST(ParameterFileMetnoKalman, description) {
      ParameterFileMetnoKalman::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
