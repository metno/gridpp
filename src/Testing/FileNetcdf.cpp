#include "../File/Arome.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

// For each test it is safe to assume that 10x10_copy.nc is identical to 10x10.nc
// After the test is done, it is safe to assume that 10x10_copy.nc is again reverted.
namespace {
   class FileNetcdf : public ::testing::Test {
      public:
         void reset10x10() const {
            Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
         };
         virtual void SetUp() {
            reset10x10();
         };
         virtual void TearDown() {
            reset10x10();
         };

      protected:
   };

   TEST_F(FileNetcdf, overwriteAttribute) {
      FileArome file = FileArome("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("history", "test512");
      EXPECT_EQ("test512", file.getGlobalAttribute("history"));
   }
   TEST_F(FileNetcdf, addAttribute) {
      FileArome file = FileArome("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("history2", "test123");
      EXPECT_EQ("test123", file.getGlobalAttribute("history2"));
   }
   TEST_F(FileNetcdf, appendAttribute) {
      // Check that appending and prepending works
      FileArome file = FileArome("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("history", "empty");
      file.prependGlobalAttribute("history",  "testing");
      file.appendGlobalAttribute("history",  "testing2");
      EXPECT_EQ("testing\nempty\ntesting2", file.getGlobalAttribute("history"));

      // Writing should not thrown an error
      std::vector<Variable::Type> vars;
      file.write(vars);
   }
   TEST_F(FileNetcdf, appendAttributeEmpty) {
      // Check that appending and prepending to an empty attribute works
      FileArome file = FileArome("testing/files/10x10_copy.nc");
      file.prependGlobalAttribute("history71623",  "value321");
      file.appendGlobalAttribute("history99311",  "value15");
      EXPECT_EQ("value321", file.getGlobalAttribute("history71623"));
      EXPECT_EQ("value15",  file.getGlobalAttribute("history99311"));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
