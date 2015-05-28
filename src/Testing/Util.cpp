#include "../Util.h"
#include "../File/Arome.h"
#include "../Variable.h"
#include "../Calibrator/Calibrator.h"
#include <algorithm>
#include <math.h>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <unistd.h>

namespace {
   class UtilTest : public ::testing::Test {
      protected:
   };

   TEST_F(UtilTest, isValid) {
      EXPECT_FALSE(Util::isValid(Util::MV));
      EXPECT_FALSE(Util::isValid(1.0/0));
      EXPECT_TRUE(Util::isValid(-15));
      EXPECT_TRUE(Util::isValid(0));
      EXPECT_TRUE(Util::isValid(3.14));
   }
   TEST_F(UtilTest, deg2rad) {
      EXPECT_FLOAT_EQ(2*Util::pi, Util::deg2rad(360));
      EXPECT_FLOAT_EQ(-2*Util::pi, Util::deg2rad(-360));
      EXPECT_FLOAT_EQ(0, Util::deg2rad(0));
      EXPECT_FLOAT_EQ(1.0/3*Util::pi, Util::deg2rad(60));
   }
   TEST_F(UtilTest, rad2deg) {
      EXPECT_FLOAT_EQ(360, Util::rad2deg(2*Util::pi));
      EXPECT_FLOAT_EQ(-360, Util::rad2deg(-2*Util::pi));
      EXPECT_FLOAT_EQ(60, Util::rad2deg(1.0/3*Util::pi));
   }
   TEST_F(UtilTest, getDistance) {
      EXPECT_FLOAT_EQ(0, Util::getDistance(60,10,60,10));
      EXPECT_FLOAT_EQ(20037508, Util::getDistance(90,10,-90,10));
      EXPECT_FLOAT_EQ(20037508, Util::getDistance(0,0,0,180));
      EXPECT_FLOAT_EQ(16879114, Util::getDistance(60.5,5.25,-84.75,-101.75));
   }
   TEST_F(UtilTest, getDate) {
      EXPECT_FLOAT_EQ(20150528, Util::getDate(1432816155));
      EXPECT_FLOAT_EQ(20001015, Util::getDate(971623971));
   }
   TEST_F(UtilTest, getDistanceInvalid) {
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(Util::MV,5.25,-84.75,-101.75));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(60.5,Util::MV,-84.75,-101.75));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(60.5,5.25,Util::MV,-101.75));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(60.5,5.25,-84.75,Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(Util::MV,5.25,Util::MV,-101.75));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(Util::MV,Util::MV,Util::MV,Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance((float) 1/0,5.25,-84.75,-101.75));
   }
   TEST_F(UtilTest, sortFirst) {
      typedef std::pair<int,int> pair;
      std::vector<pair> pairs;
      pairs.push_back(pair(0,2));
      pairs.push_back(pair(1,4));
      pairs.push_back(pair(3,5));
      pairs.push_back(pair(2,3));
      std::sort(pairs.begin(), pairs.end(), Util::sort_pair_first<int,int>());
      EXPECT_EQ(0, pairs[0].first);
      EXPECT_EQ(1, pairs[1].first);
      EXPECT_EQ(2, pairs[2].first);
      EXPECT_EQ(3, pairs[3].first);
      EXPECT_EQ(2, pairs[0].second);
      EXPECT_EQ(4, pairs[1].second);
      EXPECT_EQ(3, pairs[2].second);
      EXPECT_EQ(5, pairs[3].second);
   }
   TEST_F(UtilTest, sortSecond) {
      typedef std::pair<int,int> pair;
      std::vector<pair> pairs;
      pairs.push_back(pair(0,2));
      pairs.push_back(pair(1,4));
      pairs.push_back(pair(3,5));
      pairs.push_back(pair(2,3));
      std::sort(pairs.begin(), pairs.end(), Util::sort_pair_second<int,int>());
      EXPECT_EQ(0, pairs[0].first);
      EXPECT_EQ(2, pairs[1].first);
      EXPECT_EQ(1, pairs[2].first);
      EXPECT_EQ(3, pairs[3].first);
      EXPECT_EQ(2, pairs[0].second);
      EXPECT_EQ(3, pairs[1].second);
      EXPECT_EQ(4, pairs[2].second);
      EXPECT_EQ(5, pairs[3].second);
   }
   TEST_F(UtilTest, calcDate) {
      EXPECT_EQ(20150101, Util::calcDate(20141231, 24));
      EXPECT_EQ(20141231, Util::calcDate(20141231, 0));
      EXPECT_EQ(20141231, Util::calcDate(20141231, 1.5));
      EXPECT_EQ(20150101, Util::calcDate(20141231, 47.9));
      EXPECT_EQ(20141231, Util::calcDate(20150101, -0.1));
      EXPECT_EQ(20141231, Util::calcDate(20150101, -24));
      EXPECT_EQ(20141230, Util::calcDate(20150101, -24.01));
   }
   TEST_F(UtilTest, split) {
      std::string s = "test 1 2 3";
      std::vector<std::string> strings = Util::split(s);
      ASSERT_EQ(4, strings.size());
      EXPECT_EQ("test", strings[0]);
      EXPECT_EQ("1", strings[1]);
      EXPECT_EQ("2", strings[2]);
      EXPECT_EQ("3", strings[3]);
   }
   TEST_F(UtilTest, splitEmpty) {
      std::string s = "";
      std::vector<std::string> strings = Util::split(s);
      EXPECT_EQ(0, strings.size());
   }
   TEST_F(UtilTest, splitEmptyExtraSpace) {
      std::string s = "     ";
      std::vector<std::string> strings = Util::split(s);
      EXPECT_EQ(0, strings.size());
   }
   TEST_F(UtilTest, splitExtraSpace) {
      std::string s = "  test  1 ";
      std::vector<std::string> strings = Util::split(s);
      ASSERT_EQ(2, strings.size());
      EXPECT_EQ("test", strings[0]);
      EXPECT_EQ("1", strings[1]);
   }
   TEST_F(UtilTest, logit) {
      const float p[] = {0.2,0.9};
      const float exp[] = {-1.386294, 2.197225};
      for(int i = 0; i < sizeof(p)/sizeof(float); i++) {
         float ans = Util::logit(p[i]);
         EXPECT_FLOAT_EQ(exp[i], ans);
      }
      EXPECT_FLOAT_EQ(Util::MV, Util::logit((float) 1/0));
      EXPECT_FLOAT_EQ(Util::MV, Util::logit(-1));
      EXPECT_FLOAT_EQ(Util::MV, Util::logit(0));
      EXPECT_FLOAT_EQ(Util::MV, Util::logit(1));
   }
   TEST_F(UtilTest, invlogit) {
      const float x[] = {-2, 3};
      const float exp[] = {0.1192029, 0.9525741};
      for(int i = 0; i < sizeof(x)/sizeof(float); i++) {
         float ans = Util::invLogit(x[i]);
         EXPECT_FLOAT_EQ(exp[i], ans);
      }
      EXPECT_FLOAT_EQ(Util::MV, Util::invLogit((float) 1/0));
   }
   TEST_F(UtilTest, error) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(Util::error("test"), ".*");
   }
   TEST_F(UtilTest, warning) {
      Util::setShowWarning(false);
      Util::warning("test");
   }
   TEST_F(UtilTest, status) {
      Util::setShowStatus(false);
      Util::status("test");
   }
   TEST_F(UtilTest, clock) {
      double sleepmicros = 50000;
      double sleepsec    = sleepmicros/1e6;
      double s = Util::clock();
      usleep(sleepmicros);
      double e = Util::clock();
      EXPECT_NEAR(sleepsec, e-s, sleepsec/20);
   }
   TEST_F(UtilTest, getCurrentDate) {
      int date = Util::getCurrentDate();
      EXPECT_LT(2000000, date);
   }
   TEST_F(UtilTest, getCurrentDateString) {
      std::string date = Util::getCurrentTimeStamp();
      EXPECT_EQ(19, date.size());
   }
   TEST_F(UtilTest, copy) {
      // Check that copying a file works. Do this by modifying a file, and then overwritiing
      // it with known contents.
      bool status = Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
      EXPECT_TRUE(status);
      // Change 10x10_copy.nc
      CalibratorSmooth smooth(Variable::T);
      FileArome f1("testing/files/10x10.nc");
      FileArome f2("testing/files/10x10_copy.nc");
      FieldPtr field1 = f1.getField(Variable::T, 0);
      FieldPtr field2 = f2.getField(Variable::T, 0);
      EXPECT_EQ((*field2)(1,1,0), (*field1)(1,1,0));
      smooth.calibrate(f2);
      // Values should be different now
      EXPECT_NE((*field2)(1,1,0), (*field1)(1,1,0));
      std::vector<Variable::Type> vars(1,Variable::T);
      f2.write(vars);

      // Use copy. Values should be the same afterwards
      status = Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
      EXPECT_TRUE(status);
      FileArome f3("testing/files/10x10_copy.nc");
      FieldPtr field3 = f3.getField(Variable::T, 0);
      EXPECT_EQ(f1.getNumLat(), f2.getNumLat());
      EXPECT_EQ(f1.getNumLon(), f2.getNumLon());
      EXPECT_EQ(f1.getNumEns(), f2.getNumEns());
      EXPECT_EQ((*field3)(1,1,0), (*field1)(1,1,0));
   }
   TEST_F(UtilTest, copyInvalid) {
      bool status = Util::copy("testing/filesDoesNotExist/file1", "testing/filesDoesNotExist/file1");
      EXPECT_FALSE(status);
      status = Util::copy("testing/files/10x10.nc", "testing/filesDoesNotExist/file1");
      EXPECT_FALSE(status);
      status = Util::copy("testing/filesDoesNotExist/file1", "testing/files/10x10_copy.nc");
      EXPECT_FALSE(status);
      status = Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
      EXPECT_TRUE(status);
   }
   TEST_F(UtilTest, hasChar) {
      EXPECT_TRUE(Util::hasChar("te  2384 &$*#st", 't'));
      EXPECT_TRUE(Util::hasChar("te  2384 &$*#st", ' '));
      EXPECT_TRUE(Util::hasChar("te  2384 &$*#2t", '2'));
      EXPECT_TRUE(Util::hasChar("te  2384 &$*#2t", '&'));
      EXPECT_FALSE(Util::hasChar("te  2384 &$*#2t", ','));
      EXPECT_FALSE(Util::hasChar("te  2384 &$*#2t", 'q'));

      EXPECT_FALSE(Util::hasChar("", 'q'));
      EXPECT_FALSE(Util::hasChar("", ' '));

      EXPECT_TRUE(Util::hasChar("test\nq2", '2'));
      EXPECT_TRUE(Util::hasChar("test\nq2", 't'));
      EXPECT_TRUE(Util::hasChar("test\nq2", '\n'));
      EXPECT_FALSE(Util::hasChar("test\nq2", ' '));
      EXPECT_FALSE(Util::hasChar("test\nq2", '3'));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
