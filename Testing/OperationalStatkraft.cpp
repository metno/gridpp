#include "../File/Arome.h"
#include "../Util.h"
#include <gtest/gtest.h>

namespace {
   class OperationalStatkraftTest : public ::testing::Test {
      protected:
         OperationalStatkraftTest() {
         }
         virtual ~OperationalStatkraftTest() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
         std::string getOperationalAromeFile(int iDate, int iInit) {
            std::stringstream ss;
            ss << "/opdata/arome2_5/arome_metcoop_default2_5km_" << iDate << "_"
               << std::setfill('0') << std::setw(2) << iInit << ".nc";
            return ss.str();
         }
   };

   // Check that the AROME input files look right
   TEST_F(OperationalStatkraftTest, input) {
      int today = Util::getCurrentDate();
      int yesterday = Util::calcDate(today, -24);
      std::string filename = getOperationalAromeFile(yesterday, 0);
      FileArome file(filename, true);

      // Dimensions
      EXPECT_EQ(1, file.getNumEns());
      EXPECT_EQ(67, file.getNumTime());
      EXPECT_EQ(929, file.getNumLat());
      EXPECT_EQ(719, file.getNumLon());

      // Variables
      EXPECT_TRUE(file.hasVariable(Variable::PrecipAcc));
      EXPECT_TRUE(file.hasVariable(Variable::U));
      EXPECT_TRUE(file.hasVariable(Variable::V));
      EXPECT_TRUE(file.hasVariable(Variable::T));

      // Contents
      Field precip = *file.getField(Variable::PrecipAcc, 1);
      EXPECT_TRUE(Util::isValid(precip[0][0][0]));

   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
