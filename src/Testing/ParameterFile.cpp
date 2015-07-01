#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>
#include <boost/uuid/uuid.hpp>

namespace {
   class ParameterFileTest : public ::testing::Test {
      public:
      protected:
   };

   TEST_F(ParameterFileTest, validDownscalers) {
      ParameterFile* p0 = ParameterFile::getScheme("text", Options("file=testing/files/parameters.txt"));
      EXPECT_EQ("text", p0->name());
      EXPECT_EQ(false, p0->isLocationDependent());
      ParameterFile* p1 = ParameterFile::getScheme("text", Options("file=testing/files/parametersKriging.txt spatial=1"));
      EXPECT_EQ("text", p1->name());
      EXPECT_EQ(true, p1->isLocationDependent());
      ParameterFile* p2 = ParameterFile::getScheme("metnoKalman", Options("file=testing/files/kalmanOutput.txt"));
      EXPECT_EQ("metnoKalman", p2->name());
      ParameterFile* p3 = ParameterFile::getScheme("netcdf", Options("file=testing/files/10x10_param.nc"));
      EXPECT_EQ("netcdf", p3->name());
   }
   TEST_F(ParameterFileTest, nearestNeighbour) {
      ParameterFile* p = ParameterFile::getScheme("text", Options("file=testing/files/parametersKriging.txt spatial=1"));
      // Nearest neighbour is 5,5
      Parameters par = p->getParameters(0, Location(4.9,4.9,0));
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(4.2, par[0]);
      // Nearest neighbour at this time is 5,5
      par = p->getParameters(0, Location(-1,0.1,0));
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(4.2, par[0]);
      // Nearest neighbour at this time is 0,0
      par = p->getParameters(1, Location(-1,0.1,0));
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(4, par[0]);
   }
   TEST_F(ParameterFileTest, factoryInvalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(ParameterFile::getScheme("qwe", Options("")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
