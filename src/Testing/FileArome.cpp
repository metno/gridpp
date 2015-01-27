#include "../File/Arome.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileAromeTest : public ::testing::Test {
   };

   TEST_F(FileAromeTest, 5x5) {
      {
         FileArome from("testing/files/5x5.nc");
         FileArome to("testing/files/5x5_copy.nc");
         EXPECT_TRUE(from.hasVariable(Variable::T));
         DownscalerNearestNeighbour d(Variable::T);
         std::vector<Variable::Type> variables;
         variables.push_back(Variable::T);
         d.downscale(from, to);

         to.write(variables);
      }
      // Nearest neighbour should give the same values
      FileArome f1("testing/files/5x5.nc");
      FileArome f2("testing/files/5x5_copy.nc");
      FieldPtr p1 = f1.getField(Variable::T, 0);
      FieldPtr p2 = f2.getField(Variable::T, 0);
      EXPECT_EQ(*p1, *p2);
   }
   TEST_F(FileAromeTest, 5x5_smart) {
      {
         FileArome from("testing/files/5x5.nc");
         FileArome to("testing/files/5x5_copy.nc");
         EXPECT_TRUE(from.hasVariable(Variable::T));
         DownscalerSmart d(Variable::T);
         std::vector<Variable::Type> variables;
         variables.push_back(Variable::T);
         d.downscale(from, to);

         to.write(variables);
      }
      // Smart downscaler should give different values
      FileArome f1("testing/files/5x5.nc");
      FileArome f2("testing/files/5x5_copy.nc");
      FieldPtr p1 = f1.getField(Variable::T, 0);
      FieldPtr p2 = f2.getField(Variable::T, 0);
      EXPECT_NE(*p1, *p2);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
