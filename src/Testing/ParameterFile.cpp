#include "../ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>

namespace {
   TEST(ParameterFile, singleTime) {
      ParameterFile file("testing/files/parametersSingleTime.txt");
      ASSERT_EQ(1, file.getSize());
      Parameters par = file.getParameters(0);
      ASSERT_EQ(9, par.size());
      EXPECT_FLOAT_EQ(-1.2021, par[0]);
      EXPECT_FLOAT_EQ(0.0007985, par[8]);

      EXPECT_EQ(file.getParameters(0).getValues(), file.getParameters(10).getValues());
      EXPECT_EQ(file.getParameters(3).getValues(), file.getParameters(25).getValues());
   }

   TEST(ParameterFile, multipleTime) {
      ParameterFile file("testing/files/parametersMultipleTime.txt");
      ASSERT_EQ(8, file.getSize());
      Parameters par = file.getParameters(30);
      ASSERT_EQ(8, par.size());
      EXPECT_FLOAT_EQ(0.04198875, par[0]);
      EXPECT_FLOAT_EQ(-0.04039751, par[5]);
   }
   TEST(ParameterFile, empty) {
      ParameterFile file("testing/files/parametersEmpty.txt");
      ASSERT_EQ(6, file.getSize());
      Parameters par = file.getParameters(0);
      ASSERT_EQ(0, par.size());
   }
   TEST(ParameterFile, invalidFiles) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(ParameterFile("testing/files/parametersUnevenRows.txt"), ".*");
      EXPECT_DEATH(ParameterFile("testing/files/parametersDoesntExist.txt"), ".*");
      EXPECT_DEATH(ParameterFile("testing/files/parametersInvalidTime.txt"), ".*");
      EXPECT_DEATH(ParameterFile("testing/files/parametersInvalidEntries.txt"), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
