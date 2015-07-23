#include <gtest/gtest.h>
#include "../Variable.h"
#include "../Util.h"

namespace {

   TEST(VariableTest, invalidNames) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(Variable::getType("invalid_variable"), ".*");
      EXPECT_DEATH(Variable::getType(""), ".*");
   }
   TEST(VariableTest, conversionStringEnum) {
      std::vector<Variable::Type> variables = Variable::getAllVariables();
      for(int i = 0; i < variables.size(); i++) {
         std::string varname    = Variable::getTypeName(variables[i]);
         Variable::Type vartype = Variable::getType(varname);
         EXPECT_EQ(variables[i], vartype);
      }
   }
   TEST(VariableTest, attributes) {
      std::vector<Variable::Type> variables = Variable::getAllVariables();
      for(int i = 0; i < variables.size(); i++) {
         float min = Variable::getMin(variables[i]);
         float max = Variable::getMax(variables[i]);
         std::string units = Variable::getUnits(variables[i]);
         std::string standardNames = Variable::getStandardName(variables[i]);
      }
   }
   TEST(VariableTest, getDescriptions) {
      Variable::getDescriptions();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
