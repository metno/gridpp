#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

namespace {
   class TestDownscalerSmart : public ::testing::Test {
      public:
         vec2 makeVec2(int nLat, int nLon, const std::vector<float>& values) {
            vec2 grid;
            grid.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               grid[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  int index = i*nLon + j;
                  grid[i][j] = values[index];
               }
            }
            return grid;
         };
         // iElev must be a flat array with lon varying fastest
         void setLatLonElev(FileFake& iFile, float iLat[], float iLon[], float iElev[]) {
            vec2 lat;
            vec2 lon;
            vec2 elev;
            int nLat = iFile.getNumLat(); 
            int nLon = iFile.getNumLon();
            lat.resize(nLat);
            lon.resize(nLat);
            elev.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               lat[i].resize(nLon);
               lon[i].resize(nLon);
               elev[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  lat[i][j] = iLat[i];
                  lon[i][j] = iLon[j];
                  elev[i][j] = iElev[j+i*nLon];
               }
            }
            iFile.setLats(lat);
            iFile.setLons(lon);
            iFile.setElevs(elev);
         };
      protected:
   };

   TEST_F(TestDownscalerSmart, isValid) {
      FileFake from(3,2,1,1);
      FileFake to(1,1,1,1);
      setLatLonElev(from, (float[]) {50,55,60}, (float[]){0,10}, (float[]){3, 15, 6, 30, 20, 11});
      setLatLonElev(to,   (float[]) {54},   (float[]){9}, (float[]){10});
      int searchRadius = 10; // Search the whole grid
      int numSmart = 2;

      vec3Int I, J;
      DownscalerSmart::getSmartNeighbours(from, to, searchRadius, numSmart, I, J);

      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, I[0].size());
      ASSERT_EQ(2, I[0][0].size());

      EXPECT_EQ(2, I[0][0][0]);
      EXPECT_EQ(1, J[0][0][0]);
      EXPECT_EQ(1, I[0][0][1]);
      EXPECT_EQ(0, J[0][0][1]);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
