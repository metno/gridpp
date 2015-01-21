#include <gtest/gtest.h>
#include "../Setup.h"
#include "../Downscaler/Smart.h"
typedef Setup MetSetup;

namespace {
   class SetupTest : public ::testing::Test {
      protected:
         SetupTest() {
            // You can do set-up work for each test here.
         }

         virtual ~SetupTest() {
            // You can do clean-up work that doesn't throw exceptions here.
         }
         virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
         }

         virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
         }

   };

   TEST_F(SetupTest, test1) {
      std::vector<std::string> args;
      args.push_back("-v");
      args.push_back("T");
      args.push_back("-c");
      args.push_back("zaga");
      args.push_back("-c");
      args.push_back("accumulate");
      args.push_back("-d");
      args.push_back("smart");
      args.push_back("searchRadius");
      args.push_back("11");

      MetSetup setup;
      MetSetup::getSetup(args, setup);
      EXPECT_EQ(1, setup.variableConfigurations.size());
      EXPECT_EQ(2, setup.variableConfigurations[0].calibrators.size());
      EXPECT_EQ(Variable::T, setup.variableConfigurations[0].variable);
      EXPECT_EQ(11, ((DownscalerSmart*) setup.variableConfigurations[0].downscaler)->getSearchRadius());
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
