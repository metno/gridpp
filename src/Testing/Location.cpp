#include "../Util.h"
#include "../Location.h"
#include <gtest/gtest.h>

namespace {
   class LocationTest : public ::testing::Test {
      protected:
   };

   TEST_F(LocationTest, constructor) {
      Location loc(1,2,3);
      EXPECT_FLOAT_EQ(1, loc.lat());
      EXPECT_FLOAT_EQ(2, loc.lon());
      EXPECT_FLOAT_EQ(3, loc.elev());
   }
   TEST_F(LocationTest, setters) {
      Location loc(1,2,3);
      loc.lat(3.2);
      EXPECT_FLOAT_EQ(3.2, loc.lat());
      EXPECT_FLOAT_EQ(2, loc.lon());
      EXPECT_FLOAT_EQ(3, loc.elev());
      loc.lon(4);
      EXPECT_FLOAT_EQ(3.2, loc.lat());
      EXPECT_FLOAT_EQ(4, loc.lon());
      EXPECT_FLOAT_EQ(3, loc.elev());
      loc.elev(-3.5);
      EXPECT_FLOAT_EQ(3.2, loc.lat());
      EXPECT_FLOAT_EQ(4, loc.lon());
      EXPECT_FLOAT_EQ(-3.5, loc.elev());
   }
   TEST_F(LocationTest, order) {
      Location loc1(1,2,3);
      Location loc2(1,2,3);
      // < should be false in both cases when locations are equal
      EXPECT_FALSE(loc1 < loc2 || loc2 < loc1);
      // < should always order one high when locations are different
      loc2 = Location(1,2,4);
      EXPECT_TRUE(loc1 < loc2 || loc2 < loc1);
      loc2 = Location(1,3,2);
      EXPECT_TRUE(loc1 < loc2 || loc2 < loc1);
      loc2 = Location(2,2,3);
      EXPECT_TRUE(loc1 < loc2 || loc2 < loc1);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
