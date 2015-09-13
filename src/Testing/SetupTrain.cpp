#include <gtest/gtest.h>
#include "../SetupTrain.h"
#include "../Downscaler/Smart.h"
#include "../Calibrator/Calibrator.h"

namespace {

   TEST(TestSetupTrain, test1) {
      SetupTrain setup(Util::split("testing/files/trainingDataTemperature.txt -p text file=testing/files/calculatedParameters.txt -c gaussian -v T"));
      EXPECT_TRUE(setup.trainingData != NULL);
      EXPECT_EQ("testing/files/calculatedParameters.txt", setup.output->getFilename());
      EXPECT_EQ("gaussian", setup.method->name());
   }
   TEST(TestSetupTrain, shouldBeValid) {
      SetupTrain setup(Util::split("testing/files/trainingDataTemperature.txt -p text file=testing/files/calculatedParameters.txt -c zaga -v Precip"));
   }
   TEST(TestSetupTrain, shouldBeInValid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // No variables
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataTemperature.txt")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/10x10.nc testing/files/10x10.nc -p text file=output.txt -v")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/10x10.nc testing/files/10x10.nc -p text file=output.txt -v -c gaussian")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/10x10.nc testing/files/10x10.nc -p text file=output.txt -v -c gaussian")), ".*");
      // Invalid input file
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/weoihwoiedoiwe10x10.nc -p text file=output.txt -v T -c gaussian")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
