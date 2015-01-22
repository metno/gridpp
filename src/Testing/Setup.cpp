#include <gtest/gtest.h>
#include "../Setup.h"
#include "../Downscaler/Smart.h"
#include "../Calibrator/Calibrator.h"
typedef Setup MetSetup;

namespace {

   TEST(SetupTest, test1) {
      MetSetup setup(Util::split("input output -v T -c zaga parameters=parameters.txt -c accumulate -d smart searchRadius=11"));
      EXPECT_EQ(1,           setup.variableConfigurations.size());
      EXPECT_EQ(2,           setup.variableConfigurations[0].calibrators.size());
      EXPECT_EQ(Variable::T, setup.variableConfigurations[0].variable);
      EXPECT_EQ(11, ((DownscalerSmart*) setup.variableConfigurations[0].downscaler)->getSearchRadius());
   }
   TEST(SetupTest, variableOnly) {
      MetSetup setup(Util::split("input output -v T"));
      ASSERT_EQ(1,                          setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,                setup.variableConfigurations[0].variable);
      EXPECT_EQ(Setup::defaultDownscaler(), setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, valid) {
      MetSetup setup(Util::split("input output -v T -d smart"));
      ASSERT_EQ(1,            setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,  setup.variableConfigurations[0].variable);
      EXPECT_EQ("smart",      setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,            setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, repeatVariable) {
      MetSetup setup(Util::split("input output -v T -v T -d smart -c smooth"));
      ASSERT_EQ(1,                          setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,                setup.variableConfigurations[0].variable);
      EXPECT_EQ(Setup::defaultDownscaler(), setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, repeatDownscaler) {
      MetSetup setup(Util::split("input output -v T -d smart -d nearestNeighbour"));
      ASSERT_EQ(1,                  setup.variableConfigurations.size());
      EXPECT_EQ("nearestNeighbour", setup.variableConfigurations[0].downscaler->name());
   }
   TEST(SetupTest, complicated) {
      MetSetup setup(Util::split("input output -v T -d nearestNeighbour -d smart -c smooth -c accumulate -c smooth -v Precip -c zaga parameters=parameters.txt -d gradient"));
      ASSERT_EQ(2,            setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,  setup.variableConfigurations[0].variable);
      EXPECT_EQ("smart",      setup.variableConfigurations[0].downscaler->name());
      ASSERT_EQ(3,            setup.variableConfigurations[0].calibrators.size());
      EXPECT_EQ("smooth",     setup.variableConfigurations[0].calibrators[0]->name());
      EXPECT_EQ("accumulate", setup.variableConfigurations[0].calibrators[1]->name());
      EXPECT_EQ("smooth",     setup.variableConfigurations[0].calibrators[2]->name());

      EXPECT_EQ(Variable::Precip, setup.variableConfigurations[1].variable);
      EXPECT_EQ("gradient",   setup.variableConfigurations[1].downscaler->name());
      ASSERT_EQ(1,            setup.variableConfigurations[1].calibrators.size());
      EXPECT_EQ("zaga",       setup.variableConfigurations[1].calibrators[0]->name());
   }
   TEST(SetupTest, shouldBeValid) {
      MetSetup setup1(Util::split("input output -v T -d smart"));
      MetSetup setup2(Util::split("input output -v T -c smooth -d smart"));
      MetSetup setup3(Util::split("input output -v T -d smart -c smooth"));
      MetSetup setup4(Util::split("input output -v T -d smart -c smooth smooth"));
      MetSetup setup5(Util::split("input output -v T -d smart -c smooth smooth -v Precip -d smart"));
      MetSetup setup6(Util::split("input output -v T -d nearestNeighbour -d smart -c smooth accumulate smooth -v Precip -d smart"));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
