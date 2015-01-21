#include "../Util.h"
#include <gtest/gtest.h>

namespace {
   class UtilTest : public ::testing::Test {
      protected:
   };

   TEST_F(UtilTest, isValid) {
      EXPECT_FALSE(Util::isValid(Util::MV));
      EXPECT_FALSE(Util::isValid(1.0/0));
      EXPECT_TRUE(Util::isValid(-15));
      EXPECT_TRUE(Util::isValid(0));
      EXPECT_TRUE(Util::isValid(3.14));
   }
   TEST_F(UtilTest, deg2rad) {
      EXPECT_FLOAT_EQ(2*Util::pi, Util::deg2rad(360));
      EXPECT_FLOAT_EQ(-2*Util::pi, Util::deg2rad(-360));
      EXPECT_FLOAT_EQ(0, Util::deg2rad(0));
      EXPECT_FLOAT_EQ(1.0/3*Util::pi, Util::deg2rad(60));
   }
   TEST_F(UtilTest, rad2deg) {
      EXPECT_FLOAT_EQ(360, Util::rad2deg(2*Util::pi));
      EXPECT_FLOAT_EQ(-360, Util::rad2deg(-2*Util::pi));
      EXPECT_FLOAT_EQ(60, Util::rad2deg(1.0/3*Util::pi));
   }
   TEST_F(UtilTest, getDistance) {
      EXPECT_FLOAT_EQ(0, Util::getDistance(60,10,60,10));
      EXPECT_FLOAT_EQ(20037508, Util::getDistance(90,10,-90,10));
      EXPECT_FLOAT_EQ(20037508, Util::getDistance(0,0,0,180));
      EXPECT_FLOAT_EQ(16879114, Util::getDistance(60.5,5.25,-84.75,-101.75));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
