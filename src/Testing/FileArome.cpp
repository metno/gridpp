#include "../File/Arome.h"
#include "../Util.h"
#include <gtest/gtest.h>

namespace {
   class FileAromeTest : public ::testing::Test {
      protected:
         FileAromeTest() {
            mFilenameSmall = "testing/files/arome_small.nc";
         }
         virtual ~FileAromeTest() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
         std::string mFilenameSmall;
   };

   TEST_F(FileAromeTest, small) {
      FileArome file(mFilenameSmall);
   }
   TEST_F(FileAromeTest, operational) {
      int today = Util::getCurrentDate();
      int yesterday = Util::calcDate(today, -24);
      std::stringstream ss;
      ss << "/opdata/arome2_5/arome_metcoop_default2_5km_" << yesterday << "_00.nc";
      std::string filename = ss.str();
      std::cout << filename << std::endl;
      FileArome file(filename, true);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
