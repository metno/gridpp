#include <gtest/gtest.h>
#include "../SetupTrain.h"
#include "../Downscaler/Smart.h"
#include "../Calibrator/Calibrator.h"

namespace {

   TEST(TestSetupTrain, test1) {
      std::vector<std::string> setups;
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -c gaussian -p text file=testing/files/calculatedParameters.txt");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -p text file=testing/files/calculatedParameters.txt -c gaussian");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=testing/files/calculatedParameters.txt -v Precip -c gaussian");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=testing/files/calculatedParameters.txt -c gaussian -v Precip");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -p text file=testing/files/calculatedParameters.txt -v Precip");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -v Precip -p text file=testing/files/calculatedParameters.txt");
      for(int i = 0; i < setups.size(); i++) {
         SetupTrain setup(Util::split(setups[i]));
         ASSERT_TRUE(setup.observations.size() == 1);
         ASSERT_TRUE(setup.forecasts.size() == 1);
         EXPECT_EQ(setup.observations[0]->getFilename(), "testing/files/trainingDataObs.txt");
         EXPECT_EQ(setup.forecasts[0]->getFilename(), "testing/files/trainingDataFcst.txt");
         EXPECT_EQ("testing/files/calculatedParameters.txt", setup.output->getFilename());
         EXPECT_EQ("text", setup.output->name());
         EXPECT_EQ("gaussian", setup.method->name());
      }
   }
   TEST(TestSetupTrain, shouldBeValid) {
      // The order of -p -c -v should be irrelevant
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -c gaussian -p text file=testing/files/calculatedParameters.txt"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -p text file=testing/files/calculatedParameters.txt -c gaussian"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=testing/files/calculatedParameters.txt -v Precip -c gaussian"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=testing/files/calculatedParameters.txt -c gaussian -v Precip"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -p text file=testing/files/calculatedParameters.txt -v Precip"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -v Precip -p text file=testing/files/calculatedParameters.txt"));
   }
   TEST(TestSetupTrain, shouldBeInValid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      // Missing obs and or fcst data
      EXPECT_DEATH(SetupTrain(Util::split("")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("-v Precip -c zaga -p text file=testing/files/calculatedParameters.txt")), ".*");

      // No variables
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=output.txt")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=output.txt -v")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=output.txt -v -c gaussian")), ".*");
      // Invalid input file
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/weoihwoiedoiwe10x10.nc testing/files/oiwejiojqewqoijoijdwq.nc -p text file=output.txt -v T -c gaussian")), ".*");

      // Nothing after -p
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p")), ".*");


      // Missing -c
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=testing/files/calculatedParameters.txt -v T")), ".*");

      // Missing -p
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c zaga file=testing/files/calculatedParameters.txt -v T")), ".*");

      // Repeat arguments
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -c zaga -p text file=testing/files/calculatedParameters.txt -v T")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -v T -c zaga -p text file=testing/files/calculatedParameters.txt")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -c zaga -p text file=testing/files/calculatedParameters.txt -c gaussian")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -c zaga -c gaussian -p text file=testing/files/calculatedParameters.txt")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v Precip -c zaga -p text file=testing/files/calculatedParameters.txt -p text file=testing/files/calculatedParameters.txt")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text file=testing/files/calculatedParameters.txt -v Precip -c zaga -p text file=testing/files/calculatedParameters.txt")), ".*");

      // Junk after variable
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p text -v Precip file=testing/files/calculatedParameters.txt -c zaga")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
