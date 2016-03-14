#include "../KDTree.h"
#include "../File/Arome.h"
#include "../Variable.h"
#include "../Calibrator/Calibrator.h"
#include <algorithm>
#include <math.h>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <unistd.h>

namespace {
   class KDTreeTest : public ::testing::Test {
      protected:
   };
   // No locations in tree
   TEST_F(KDTreeTest, empty) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      vec2 lats, lons;
      EXPECT_DEATH(KDTree(lats, lons), ".*");
   }
   // One point at 3,2
   TEST_F(KDTreeTest, single) {
      vec2 lats, lons;
      std::vector<float> lat(1,3), lon(1,2);
      lats.push_back(lat);
      lons.push_back(lon);
      KDTree tree(lats, lons);
      int I, J;
      tree.getNearestNeighbour(3,2, I, J);
      EXPECT_EQ(0, I);
      EXPECT_EQ(0, J);
      tree.getNearestNeighbour(2,1, I, J);
      EXPECT_EQ(0, I);
      EXPECT_EQ(0, J);
   }
   // A single row of lat/lon locations: [0,3] [2,0], [0,0], [2,2]
   TEST_F(KDTreeTest, 1row) {
      vec2 lats, lons;
      std::vector<float> lat(4,0), lon(4,0);
      lat[0] = 3; lon[0] = 3;
      lat[1] = 2; lon[1] = 0;
      lat[2] = 0; lon[2] = 0;
      lat[3] = 2; lon[3] = 2;
      lats.push_back(lat);
      lons.push_back(lon);
      KDTree tree(lats, lons);
      int I, J;
      tree.getNearestNeighbour(3,3, I, J);
      EXPECT_EQ(0, I);
      EXPECT_EQ(0, J);
      tree.getNearestNeighbour(0.5,0.9, I, J);
      EXPECT_EQ(0, I);
      EXPECT_EQ(2, J);
      // Either point 0,1 or 0,3 is the nearest
      tree.getNearestNeighbour(2.1,-0.1, I, J);
      EXPECT_EQ(0, I);
      EXPECT_TRUE(J == 1 || J == 3);
   }
   // Two rows of lat/lon locations:
   // [0,0] [0,1] [0,2] [0,3]
   // [1,0] [1,1] [1,2] [1,3]
   TEST_F(KDTreeTest, matrix) {
      vec2 lats, lons;
      std::vector<float> lat(4,0), lon(4,0);
      lat[0] = 0; lon[0] = 0;
      lat[1] = 0; lon[1] = 1;
      lat[2] = 0; lon[2] = 2;
      lat[3] = 0; lon[3] = 3;
      lats.push_back(lat);
      lons.push_back(lon);
      lat[0] = 1; lon[0] = 0;
      lat[1] = 1; lon[1] = 1;
      lat[2] = 1; lon[2] = 2;
      lat[3] = 1; lon[3] = 3;
      lats.push_back(lat);
      lons.push_back(lon);
      KDTree tree(lats, lons);
      int I, J;
      tree.getNearestNeighbour(0,0, I, J);
      EXPECT_EQ(0, I);
      EXPECT_EQ(0, J);
      tree.getNearestNeighbour(1.1, 0.6, I, J);
      EXPECT_EQ(1, I);
      EXPECT_EQ(1, J);
      tree.getNearestNeighbour(0.2, 2.4, I, J);
      EXPECT_EQ(0, I);
      EXPECT_EQ(2, J);
      tree.getNearestNeighbour(10, 10, I, J);
      EXPECT_EQ(1, I);
      EXPECT_EQ(3, J);
      tree.getNearestNeighbour(-10, 10, I, J);
      EXPECT_EQ(0, I);
      EXPECT_EQ(3, J);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
