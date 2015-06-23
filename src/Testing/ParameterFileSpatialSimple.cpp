#include "../ParameterFile/ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   TEST(ParameterFileSpatialSimple, default) {
      ParameterFileSpatialSimple file("testing/files/parametersKriging.txt");
      std::vector<Location> locations = file.getLocations();
      ASSERT_EQ(4, locations.size());
      std::set<Location> locationsSet(locations.begin(), locations.end());
      std::set<Location> expected;
      expected.insert(Location(9, 9, 800));
      expected.insert(Location(9, 8, 850));
      expected.insert(Location(5, 5, 140));
      expected.insert(Location(0, 0, 150));

      EXPECT_EQ(locationsSet, expected);
   }
   TEST(ParameterFileSpatialSimple, write) {
      ParameterFileSpatialSimple file("testing/files/parametersKriging.txt");
      file.write("testing/files/parametersKriging2.txt");
      ParameterFileSpatialSimple file2("testing/files/parametersKriging2.txt");
      // Should have the same locations
      std::vector<Location> locations1 = file.getLocations();
      std::vector<Location> locations2 = file.getLocations();
      std::sort(locations1.begin(), locations1.end());
      std::sort(locations2.begin(), locations2.end());
      EXPECT_EQ(locations1, locations2);
   }
   TEST(ParameterFileSpatialSimple, emptyFile) {
      ParameterFileSpatialSimple file("testing/files/89h9382he9823he92.txt");
      std::vector<Location> locations = file.getLocations();
      EXPECT_EQ(0, locations.size());
   }
   TEST(ParameterFileSpatialSimple, writeError) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      ParameterFileSpatialSimple file("testing/files/parametersKriging.txt");
      // Shouldn't be able to write to a directory
      EXPECT_DEATH(file.write("testing/files/"), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
