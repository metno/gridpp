#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Smart.h"
#include <gtest/gtest.h>

namespace {
   class TestDownscalerSmart : public ::testing::Test {
      protected:
   };

   TEST_F(TestDownscalerSmart, isValid) {
      FileFake from(2,2,1,1);
      FileFake to(2,2,1,1);
      vec2 fromLat(2);
      fromLat[0].resize(2);
      fromLat[1].resize(2);
      fromLat[0][0] = 0
      vec3Int I;
      vec3Int J;
      DownscalerSmart::getSmartNeighbours(from, to, 10, 5, I, J);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
