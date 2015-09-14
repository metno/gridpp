#include "../File/Arome.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileAromeTest : public ::testing::Test {
   };

   TEST_F(FileAromeTest, 10x10) {
      {
         FileArome from("testing/files/10x10.nc");
         FileArome to("testing/files/10x10_copy.nc");
         EXPECT_TRUE(from.hasVariable(Variable::T));
         EXPECT_FALSE(from.hasVariable(Variable::Pop6h));
         DownscalerNearestNeighbour d(Variable::T, Options());
         std::vector<Variable::Type> variables;
         variables.push_back(Variable::T);
         d.downscale(from, to);

         to.write(variables);
      }
      // Nearest neighbour should give the same values
      FileArome f1("testing/files/10x10.nc");
      FileArome f2("testing/files/10x10_copy.nc");
      FieldPtr p1 = f1.getField(Variable::T, 0);
      FieldPtr p2 = f2.getField(Variable::T, 0);
      EXPECT_EQ(*p1, *p2);
   }
   TEST_F(FileAromeTest, 10x10_smart) {
      // Smart downscaler should give different values
      {
         FileArome from("testing/files/10x10.nc");
         FileArome to("testing/files/10x10_copy.nc");
         EXPECT_TRUE(from.hasVariable(Variable::T));
         DownscalerSmart d(Variable::T, Options());
         std::vector<Variable::Type> variables;
         variables.push_back(Variable::T);
         d.downscale(from, to);

         to.write(variables);
      }
      FileArome f1("testing/files/10x10.nc");
      FileArome f2("testing/files/10x10_copy.nc");
      FieldPtr p1 = f1.getField(Variable::T, 0);
      FieldPtr p2 = f2.getField(Variable::T, 0);
      EXPECT_NE(*p1, *p2);
   }
   TEST_F(FileAromeTest, variables) {
      FileArome file("testing/files/10x10.nc");
      std::vector<Variable::Type> variables = Variable::getAllVariables();
      for(int i = 0; i < variables.size(); i++) {
         std::string var = file.getVariableName(variables[i]);
      }
   }
   TEST_F(FileAromeTest, validFiles) {
      FileArome file1("testing/files/validArome1.nc");
      FileArome file2("testing/files/validArome2.nc");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
