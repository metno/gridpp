#include <gtest/gtest.h>
#include "../Setup.h"
#include "../Downscaler/Smart.h"
#include "../Calibrator/Calibrator.h"
typedef Setup MetSetup;

namespace {

   TEST(SetupTest, test1) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -c zaga parameters=testing/files/parameters.txt -c accumulate -d smart searchRadius=11"));
      EXPECT_EQ(1,           setup.variableConfigurations.size());
      EXPECT_EQ(2,           setup.variableConfigurations[0].calibrators.size());
      EXPECT_EQ(Variable::T, setup.variableConfigurations[0].variable);
      EXPECT_EQ(11, ((DownscalerSmart*) setup.variableConfigurations[0].downscaler)->getSearchRadius());
   }
   TEST(SetupTest, variableOnly) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T"));
      ASSERT_EQ(1,                          setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,                setup.variableConfigurations[0].variable);
      EXPECT_EQ(Setup::defaultDownscaler(), setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, valid) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart"));
      ASSERT_EQ(1,            setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,  setup.variableConfigurations[0].variable);
      EXPECT_EQ("smart",      setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,            setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, repeatVariable) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -v T -d smart -c neighbourhood"));
      ASSERT_EQ(1,                          setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,                setup.variableConfigurations[0].variable);
      EXPECT_EQ(Setup::defaultDownscaler(), setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, repeatDownscaler) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart -d nearestNeighbour"));
      ASSERT_EQ(1,                  setup.variableConfigurations.size());
      EXPECT_EQ("nearestNeighbour", setup.variableConfigurations[0].downscaler->name());
   }
   TEST(SetupTest, complicated) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d nearestNeighbour -d smart -c neighbourhood -c accumulate -c neighbourhood -v Precip -c zaga parameters=testing/files/parameters.txt -d gradient"));
      ASSERT_EQ(2,            setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];
      EXPECT_EQ(Variable::T,  varconf.variable);
      EXPECT_EQ("smart",      varconf.downscaler->name());
      ASSERT_EQ(3,            varconf.calibrators.size());
      EXPECT_EQ("neighbourhood",     varconf.calibrators[0]->name());
      EXPECT_EQ("accumulate", varconf.calibrators[1]->name());
      EXPECT_EQ("neighbourhood",     varconf.calibrators[2]->name());

      EXPECT_EQ(Variable::Precip, setup.variableConfigurations[1].variable);
      EXPECT_EQ("gradient",   setup.variableConfigurations[1].downscaler->name());
      ASSERT_EQ(1,            setup.variableConfigurations[1].calibrators.size());
      EXPECT_EQ("zaga",       setup.variableConfigurations[1].calibrators[0]->name());
   }
   TEST(SetupTest, variableOptionsSingle) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T write=0"));
      ASSERT_EQ(1,            setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];

      Options vOptions = varconf.variableOptions;
      bool doWrite = true;
      ASSERT_TRUE(vOptions.getValue("write", doWrite));
      ASSERT_FALSE(vOptions.getValue("-d", doWrite));
      EXPECT_EQ(0, doWrite);
   }
   TEST(SetupTest, variableOptionsMultiple) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -v P write=0 -v RH -v U test=2 -d smart -v V -v Precip new=2.1 -c neighbourhood"));
      ASSERT_EQ(6,            setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];

      EXPECT_EQ(Variable::T, varconf.variable);
      Options vOptions = varconf.variableOptions;
      bool doWrite = true;
      float value = -1;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[1];
      EXPECT_EQ(Variable::P, varconf.variable);
      vOptions = varconf.variableOptions;
      ASSERT_TRUE(vOptions.getValue("write", doWrite));
      EXPECT_EQ(0, doWrite);
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[2];
      EXPECT_EQ(Variable::RH, varconf.variable);
      vOptions = varconf.variableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[3];
      EXPECT_EQ(Variable::U, varconf.variable);
      EXPECT_EQ("smart", varconf.downscaler->name());
      vOptions = varconf.variableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      ASSERT_TRUE(vOptions.getValue("test", value));
      EXPECT_FLOAT_EQ(2, value);
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[4];
      EXPECT_EQ(Variable::V, varconf.variable);
      vOptions = varconf.variableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[5];
      EXPECT_EQ(Variable::Precip, varconf.variable);
      vOptions = varconf.variableOptions;
      ASSERT_TRUE(vOptions.getValue("new", value));
      EXPECT_FLOAT_EQ(2.1, value);
      EXPECT_FALSE(vOptions.getValue("-v", value));
      EXPECT_EQ(1, varconf.calibrators.size());
      EXPECT_EQ("neighbourhood", varconf.calibrators[0]->name());
   }
   TEST(SetupTest, shouldBeValid) {
      MetSetup setup1(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart"));
      MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -c neighbourhood -d smart"));
      MetSetup setup3(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart -c neighbourhood"));
      MetSetup setup4(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart -c neighbourhood neighbourhood"));
      MetSetup setup5(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart -c neighbourhood neighbourhood -v Precip -d smart"));
      MetSetup setup6(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d nearestNeighbour -d smart -c neighbourhood accumulate neighbourhood -v Precip -d smart"));
      MetSetup setup7(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d nearestNeighbour -v Precip -d smart"));
      MetSetup setup8(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart numSmart=2 -c neighbourhood -v Precip -d smart"));
      MetSetup setup9(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart numSmart=2 -v Precip -d smart"));
      MetSetup setup10(Util::split("testing/files/10x10.nc testing/files/10x10_copy.nc -v T -d smart numSmart=2 -v Precip -d smart"));
   }
   TEST(SetupTest, shouldBeInValid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // No variables
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc")), ".*");
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v")), ".*");
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v -d smart")), ".*");
      // Too many files (should have max two files)
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v -d smart testing/files/10x10.nc")), ".*");
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc testing/files/10x10.nc -v -d smart testing/files/10x10.nc")), ".*");
      // Invalid input file
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/weoihwoiedoiwe10x10.nc testing/files/10x10.nc -v -d smart testing/files/10x10.nc")), ".*");
      // Invalid output file
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/weoihwoiedoiwe10x10.nc -v -d smart testing/files/10x10.nc")), ".*");
      // Nothing after downscaler
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v Precip -d")), ".*");
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v Precip -d -c neighbourhood")), ".*");
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v Precip -c neighbourhood -d")), ".*");
      // Nothing after calibrator
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v Precip -c")), ".*");
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v Precip -d nearest -c")), ".*");
      EXPECT_DEATH(MetSetup setup2(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v Precip -c -d nearest")), ".*");
   }
   TEST(SetupTest, defaultDownscaler) {
      std::string downscaler = Setup::defaultDownscaler();
      EXPECT_NE("", downscaler);
   }
   TEST(SetupTest, variableConfiguration) {
      VariableConfiguration varconf;
   }
   TEST(SetupTest, destructor) {
      {
         MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v T -d smart numSmart=2 -v Precip -d smart"));
      }
   }
   TEST(SetupTest, inputoutputOptions) {
      MetSetup setup0(Util::split("testing/files/10x10.nc option1=1 testing/files/10x10.nc option2=2 -v T write=1 -d smart numSmart=2"));
      int i;
      EXPECT_TRUE(setup0.inputOptions.getValue("option1", i));
      EXPECT_FALSE(setup0.inputOptions.getValue("option2", i));
      EXPECT_FALSE(setup0.inputOptions.getValue("write", i));
      EXPECT_EQ(1, i);
      EXPECT_TRUE(setup0.outputOptions.getValue("option2", i));
      EXPECT_FALSE(setup0.outputOptions.getValue("option1", i));
      EXPECT_FALSE(setup0.outputOptions.getValue("write", i));
      EXPECT_EQ(2, i);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
