#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

namespace {
   class TestDownscaler : public ::testing::Test {
      public:
      protected:
   };

   TEST_F(TestDownscaler, valid) {
      Downscaler* d0 = Downscaler::getScheme("nearestNeighbour", Variable::T, Options());
      Downscaler* d1 = Downscaler::getScheme("smart", Variable::T, Options("searchRadius=3 numSmart=2"));
      Downscaler* d2 = Downscaler::getScheme("gradient", Variable::T, Options("searchRadius=5 constantGradient=0.04"));
      EXPECT_EQ(3, ((DownscalerSmart*) d1)->getSearchRadius());
      EXPECT_EQ(2, ((DownscalerSmart*) d1)->getNumSmart());
      EXPECT_EQ(5, ((DownscalerGradient*) d2)->getSearchRadius());
      EXPECT_FLOAT_EQ(0.04, ((DownscalerGradient*) d2)->getConstantGradient());
   }
   TEST_F(TestDownscaler, inalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(Downscaler::getScheme("woehiowciwofew", Variable::T, Options("searchRadius=5")), ".*");
      EXPECT_DEATH(Downscaler::getScheme("woehiowciwofew", Variable::T, Options()), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
