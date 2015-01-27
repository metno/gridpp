#include "../Parameters.h"
#include <algorithm>
#include <math.h>
#include <gtest/gtest.h>
#include "../Util.h"

namespace {
   class ParametersTest : public ::testing::Test {
      protected:
   };

   TEST_F(ParametersTest, empty) {
      Parameters parameters;
      EXPECT_EQ(0, parameters.size());
      std::vector<float> values = parameters.getValues();
      EXPECT_EQ(0, values.size());
   }
   TEST_F(ParametersTest, empty2) {
      std::vector<float> values;
      Parameters parameters;
      EXPECT_EQ(0, parameters.size());
      std::vector<float> values2 = parameters.getValues();
      EXPECT_EQ(0, values2.size());
   }
   TEST_F(ParametersTest, access) {
      std::vector<float> values(3);
      values[0] = 2;
      values[1] = 3.3;
      values[2] = 0;
      Parameters par(values);
      EXPECT_EQ(3, par.size());
      EXPECT_FLOAT_EQ(2, par[0]);
      EXPECT_FLOAT_EQ(3.3, par[1]);
      EXPECT_FLOAT_EQ(0, par[2]);

      std::vector<float> values2 = par.getValues();
      EXPECT_EQ(3, values2.size());
      EXPECT_FLOAT_EQ(2, values2[0]);
      EXPECT_FLOAT_EQ(3.3, values2[1]);
      EXPECT_FLOAT_EQ(0, values2[2]);
   }
   TEST_F(ParametersTest, emptyAccess) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      Parameters par;
      EXPECT_DEATH(par[-1], ".*");
      EXPECT_DEATH(par[0], ".*");
      EXPECT_DEATH(par[1], ".*");
   }
   TEST_F(ParametersTest, invalidAccess) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      std::vector<float> values(3);
      values[0] = 2;
      values[1] = 3.3;
      values[2] = 0;
      Parameters par(values);

      EXPECT_DEATH(par[-1], ".*");
      EXPECT_DEATH(par[Util::MV], ".*");
      EXPECT_DEATH(par[3], ".*");
      EXPECT_DEATH(par[100], ".*");
   }
   TEST_F(ParametersTest, constAccess) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      std::vector<float> values(3);
      values[0] = 2;
      values[1] = 3.3;
      values[2] = 0;
      Parameters par(values);
      const Parameters par2 = par;
      EXPECT_FLOAT_EQ(2, par2[0]);
      EXPECT_FLOAT_EQ(3.3, par2[1]);
      EXPECT_FLOAT_EQ(0, par2[2]);
      EXPECT_DEATH(par2[3], ".*");
      EXPECT_DEATH(par2[-1], ".*");
   }
   TEST_F(ParametersTest, emptyConstAccess) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      Parameters par;
      const Parameters par2 = par;
      EXPECT_DEATH(par2[-1], ".*");
      EXPECT_DEATH(par2[0], ".*");
      EXPECT_DEATH(par2[1], ".*");
   }
   TEST_F(ParametersTest, assignment) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      std::vector<float> values(3);
      values[0] = 2;
      values[1] = 3.3;
      values[2] = 0;
      Parameters par(values);
      par[0] = 4;
      par[1] = 1;
      EXPECT_FLOAT_EQ(4, par[0]);
      EXPECT_FLOAT_EQ(1, par[1]);
      EXPECT_FLOAT_EQ(0, par[2]);
      EXPECT_DEATH(par[3], ".*");
      EXPECT_DEATH(par[-1], ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
