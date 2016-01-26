#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

namespace {
   class ParameterFileTest : public ::testing::Test {
      public:
         Parameters createParameters(float i1, float i2, float i3) {
            std::vector<float> values;
            values.push_back(i1);
            values.push_back(i2);
            values.push_back(i3);
            return Parameters(values);
         };
      protected:
   };

   TEST_F(ParameterFileTest, validDownscalers) {
      ParameterFile* p0 = ParameterFile::getScheme("text", Options("file=testing/files/parameters.txt"));
      EXPECT_EQ("text", p0->name());
      EXPECT_FALSE(p0->isLocationDependent());
      ParameterFile* p1 = ParameterFile::getScheme("text", Options("file=testing/files/parametersKriging.txt spatial=1"));
      EXPECT_EQ("text", p1->name());
      EXPECT_TRUE(p1->isLocationDependent());
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

      // Nearest neighbour for 9,9 is 9,9 for time 0 and 5,5 for time 1
      Location loc(Util::MV, Util::MV);
      p->getNearestLocation(0, Location(9,9,0), loc);
      EXPECT_FLOAT_EQ(9, loc.lat());
      EXPECT_FLOAT_EQ(9, loc.lon());
      p->getNearestLocation(1, Location(9,9,0), loc);
      EXPECT_FLOAT_EQ(5, loc.lat());
      EXPECT_FLOAT_EQ(5, loc.lon());

      par = p->getParameters(0, Location(9,9,0));
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(-1.1, par[0]);
      par = p->getParameters(1, Location(9,9,0));
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(-5.4, par[0]);
   }
   TEST_F(ParameterFileTest, factoryInvalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(ParameterFile::getScheme("qwe", Options("")), ".*");
   }
   TEST_F(ParameterFileTest, setParameters) {
      ParameterFile* p = ParameterFile::getScheme("text", Options("file=testing/files/temp1231.txt spatial=1"));
      Parameters par = createParameters(5,2,3);
      Location loc = Location(0,0,0);
      Parameters par2;
      p->setParameters(createParameters(1,2,3),    0, Location(0,0,0));
      p->setParameters(createParameters(4,5,6),    3, Location(0,0,0));
      p->setParameters(createParameters(7,8,9),    2, Location(1,0,0));
      p->setParameters(createParameters(10,11,12), 3, Location(1,0,0));
      // Location 0,0,0 for time 0
      par2 = p->getParameters(0, loc);
      ASSERT_EQ(3, par2.size());
      EXPECT_FLOAT_EQ(1, par2[0]);
      // No parametesr for time 1
      par2 = p->getParameters(1, loc);
      ASSERT_EQ(0, par2.size());
      // Location 1,0,0 for time 2
      par2 = p->getParameters(2, loc);
      ASSERT_EQ(3, par2.size());
      EXPECT_FLOAT_EQ(7, par2[0]);
      // Location 1,0,0 for time 2
      par2 = p->getParameters(2, loc, false);
      ASSERT_EQ(0, par2.size());
      // Location 0,0,0 for time 3
      par2 = p->getParameters(3, loc);
      ASSERT_EQ(3, par2.size());
      EXPECT_FLOAT_EQ(4, par2[0]);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
