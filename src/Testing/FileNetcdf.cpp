#include "../File/Arome.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileNetcdf : public ::testing::Test {
      public:
         void reset10x10() const {
            Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
         };
      protected:
   };

   TEST_F(FileNetcdf, overwriteAttribute) {
      reset10x10();
      FileArome file = FileArome("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("history",  "testing");
      std::vector<Variable::Type> vars;
      file.write(vars);
   }
   TEST_F(FileNetcdf, addAttribute) {
      reset10x10();
      FileArome file = FileArome("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("history2",  "testing");
      std::vector<Variable::Type> vars;
      file.write(vars);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
