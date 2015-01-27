#include "../File/Arome.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileTest : public ::testing::Test {
   };

   TEST_F(FileTest, 5x5) {
      File* file = File::getScheme("testing/files/5x5.nc");
      EXPECT_EQ("arome", ((FileArome*)file)->name());
   }
   TEST_F(FileTest, 5x5_smart) {
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
      // Nearest neighbour should give the same values
      FileArome f1("testing/files/5x5.nc");
      FileArome f2("testing/files/5x5_copy.nc");
      FieldPtr p1 = f1.getField(Variable::T, 0);
      FieldPtr p2 = f2.getField(Variable::T, 0);
      EXPECT_NE(*p1, *p2);
   }
   TEST_F(FileTest, deriveVariables) {
      FileArome from("testing/files/5x5.nc");
      EXPECT_TRUE(from.hasVariable(Variable::Precip));
      FieldPtr precip = from.getField(Variable::Precip, 0);
      FieldPtr precipAcc = from.getField(Variable::PrecipAcc, 0);
      EXPECT_FLOAT_EQ(0.911191, (*precip)[5][5][0]);
      EXPECT_FLOAT_EQ(0,        (*precipAcc)[5][5][0]);
   }
   TEST_F(FileTest, hasSameDimensions) {
      FileArome f1("testing/files/5x5.nc");
      FileArome f2("testing/files/5x5_copy.nc");
      FileFake f3(3,3,1,1);
      EXPECT_TRUE(f1.hasSameDimensions(f2));
      EXPECT_TRUE(f2.hasSameDimensions(f1));
      EXPECT_FALSE(f1.hasSameDimensions(f3));
      EXPECT_FALSE(f3.hasSameDimensions(f1));
   }
   TEST_F(FileTest, initNewVariable) {
      FileArome f1("testing/files/5x5.nc");
      EXPECT_FALSE(f1.hasVariable(Variable::U));
      f1.initNewVariable(Variable::U);
      EXPECT_TRUE(f1.hasVariable(Variable::U));
      FieldPtr field = f1.getField(Variable::U, 0);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
