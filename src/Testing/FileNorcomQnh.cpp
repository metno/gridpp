#include "../File/NorcomQnh.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileNorcomQnhTest : public ::testing::Test {
   };

   TEST_F(FileNorcomQnhTest, options) {
      FileNorcomQnh file("testing/files/test.txt", Options("lats=1,2 lons=2,3 elevs=100,120 names=point1,point2 numTimes=2 startTime=0 endTime=2"));
      vec2 lats = file.getLats();
      vec2 lons = file.getLons();
      vec2 elevs = file.getElevs();
      ASSERT_EQ(1, lats.size());
      ASSERT_EQ(1, lons.size());
      ASSERT_EQ(1, elevs.size());
      ASSERT_EQ(2, lats[0].size());
      ASSERT_EQ(2, lons[0].size());
      ASSERT_EQ(2, elevs[0].size());
      EXPECT_FLOAT_EQ(1, lats[0][0]);
      EXPECT_FLOAT_EQ(2, lats[0][1]);
      EXPECT_FLOAT_EQ(2, lons[0][0]);
      EXPECT_FLOAT_EQ(3, lons[0][1]);
      EXPECT_FLOAT_EQ(100, elevs[0][0]);
      EXPECT_FLOAT_EQ(120, elevs[0][1]);
   }
   TEST_F(FileNorcomQnhTest, asOutput) {
      FileArome from("testing/files/10x10.nc");
      FileNorcomQnh to("testing/files/test.txt", Options("lats=1,2 lons=2,3 elevs=100,120 names=point1,point2 numTimes=2 startTime=0 endTime=1"));
      DownscalerNearestNeighbour d = DownscalerNearestNeighbour(Variable::P);
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      std::vector<Variable::Type> variables(1, Variable::P);
      to.write(variables);

      // FieldPtr field = to.getField(Variable::P, 0);
      // EXPECT_FLOAT_EQ(1, (*field)(0,0,0));
   }
   TEST_F(FileNorcomQnhTest, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // lats, lons, elevs not the same
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1,2 lons=2 elevs=3 numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=2 lons=2,3,2 elevs=3 numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=2 lons=2 elevs=3,2 numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=2 elevs=3 names=q,w numTimes=2 startTime=0 endTime=1")), ".*");

      // Start/end time not in ascending order
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=2 elevs=3 names=q numTimes=2 startTime=1 endTime=0")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
