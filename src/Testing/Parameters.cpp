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
      std::vector<float> values;
      values.push_back(3.154);
      values.push_back(2);
      Parameters parameters(values);
      EXPECT_EQ(2, parameters.size());
      EXPECT_FLOAT_EQ(3.154, parameters[0]);
      EXPECT_FLOAT_EQ(2, parameters[1]);

      std::vector<float> values2 = parameters.getValues();
      EXPECT_EQ(2, values2.size());
      EXPECT_FLOAT_EQ(3.154, values2[0]);
      EXPECT_FLOAT_EQ(2, values2[1]);
   }
   TEST_F(ParametersTest, access2) {
      std::vector<float> values(3);
      values[0] = 2;
      values[1] = 3.3;
      values[2] = 0;
      Parameters par(values);
      EXPECT_FLOAT_EQ(2, par[0]);
      EXPECT_FLOAT_EQ(3.3, par[1]);
      EXPECT_FLOAT_EQ(0, par[2]);
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
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
