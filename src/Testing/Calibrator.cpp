#include "../Calibrator/Calibrator.h"
#include "../Util.h"
#include "../Options.h"
#include <gtest/gtest.h>
#include <vector>

namespace {
   class TestCalibrator : public ::testing::Test {
      public:
         std::vector<float> vec(float iVector[], int iSize) {
            std::vector<float> vec(iSize, 0);
            for(int i = 0; i < iSize; i++) {
               vec[i] = iVector[i];
            }
            return vec;
         };
      protected:
   };

   TEST_F(TestCalibrator, shuffle) {
      const std::vector<float> vec1 = vec((float[]) {5,1,4,7,6,2,3}, 7);
            std::vector<float> vec2 = vec((float[]) {32,14,21,0,11,2,5}, 7);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(7, vec1.size());
      ASSERT_EQ(7, vec2.size());
      EXPECT_EQ(14, vec2[0]);
      EXPECT_EQ(0, vec2[1]);
      EXPECT_EQ(11, vec2[2]);
      EXPECT_EQ(32, vec2[3]);
      EXPECT_EQ(21, vec2[4]);
      EXPECT_EQ(2, vec2[5]);
      EXPECT_EQ(5, vec2[6]);
   }
   TEST_F(TestCalibrator, shuffleDifferentSizes) {
      const std::vector<float> vec1 = vec((float[]) {5,4}, 2);
            std::vector<float> vec2 = vec((float[]) {1,2,3}, 3);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(2, vec1.size());
      ASSERT_EQ(3, vec2.size());
      // Should remain unchanged
      EXPECT_EQ(1, vec2[0]);
      EXPECT_EQ(2, vec2[1]);
      EXPECT_EQ(3, vec2[2]);
   }
   TEST_F(TestCalibrator, shuffleZeroSize) {
      const std::vector<float> vec1;
            std::vector<float> vec2 = vec((float[]) {1,2,3}, 3);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(0, vec1.size());
      ASSERT_EQ(3, vec2.size());
      // Should remain unchanged
      EXPECT_EQ(1, vec2[0]);
      EXPECT_EQ(2, vec2[1]);
      EXPECT_EQ(3, vec2[2]);
   }
   TEST_F(TestCalibrator, shuffleZeroSize2) {
      const std::vector<float> vec1;
            std::vector<float> vec2;
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(0, vec1.size());
      ASSERT_EQ(0, vec2.size());
   }
   TEST_F(TestCalibrator, shuffleWithMissing) {
      const std::vector<float> vec1 = vec((float[]) {3, Util::MV, 19,3}, 4);
            std::vector<float> vec2 = vec((float[]) {1,2,4,3}, 4);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(4, vec1.size());
      ASSERT_EQ(4, vec2.size());
      // Should remain unchanged
      EXPECT_EQ(1, vec2[0]);
      EXPECT_EQ(2, vec2[1]);
      EXPECT_EQ(4, vec2[2]);
      EXPECT_EQ(3, vec2[3]);
   }
   TEST_F(TestCalibrator, shuffleWithMissing2) {
      const std::vector<float> vec1 = vec((float[]) {3, 1, 19,3}, 4);
            std::vector<float> vec2 = vec((float[]) {1,2,Util::MV,3}, 4);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(4, vec1.size());
      ASSERT_EQ(4, vec2.size());
      // Should remain unchanged
      EXPECT_EQ(1, vec2[0]);
      EXPECT_EQ(2, vec2[1]);
      EXPECT_EQ(Util::MV, vec2[2]);
      EXPECT_EQ(3, vec2[3]);
   }
   TEST_F(TestCalibrator, shuffleRepeated) {
      const std::vector<float> vec1 = vec((float[]) {3,1,7,1}, 4);
            std::vector<float> vec2 = vec((float[]) {1,2,4,3}, 4);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(4, vec1.size());
      ASSERT_EQ(4, vec2.size());
      EXPECT_EQ(3, vec2[0]);
      EXPECT_TRUE(vec2[1]==1 || vec2[1]==2);
      EXPECT_EQ(4, vec2[2]);
      EXPECT_TRUE(vec2[3]==1 || vec2[3]==2);
   }
   TEST_F(TestCalibrator, shuffleRepeated2) {
      const std::vector<float> vec1 = vec((float[]) {3,2,7,1}, 4);
            std::vector<float> vec2 = vec((float[]) {1,2,1,3}, 4);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(4, vec1.size());
      ASSERT_EQ(4, vec2.size());
      EXPECT_EQ(2, vec2[0]);
      EXPECT_TRUE(vec2[1]==1 || vec2[1]==2);
      EXPECT_EQ(3, vec2[2]);
      EXPECT_TRUE(vec2[3]==1 || vec2[3]==2);
   }
   TEST_F(TestCalibrator, factoryZaga) {
      Calibrator* c;
      c = Calibrator::getScheme("zaga", Options("variable=T parameters=testing/files/parameters.txt fracThreshold=0.9"));
      EXPECT_TRUE(c);
      EXPECT_EQ("zaga", c->name());
      EXPECT_FLOAT_EQ(0.9, ((CalibratorZaga*) c)->getFracThreshold());
      delete c;
   }
   TEST_F(TestCalibrator, factoryNeighbourhood) {
      Calibrator* c;
      c = Calibrator::getScheme("neighbourhood", Options("variable=Precip radius=3"));
      EXPECT_TRUE(c);
      EXPECT_EQ("neighbourhood", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryPhase) {
      Calibrator* c;
      c = Calibrator::getScheme("phase", Options("variable=Precip parameters=testing/files/parametersPhase.txt minPrecip=0.771 useWetbulb=0"));
      EXPECT_TRUE(c);
      EXPECT_EQ("phase", c->name());
      EXPECT_FLOAT_EQ(0.771, ((CalibratorPhase*) c)->getMinPrecip());
      EXPECT_FALSE(((CalibratorPhase*) c)->getUseWetbulb());
      delete c;
   }
   TEST_F(TestCalibrator, factoryRegression) {
      Calibrator* c;
      c = Calibrator::getScheme("regression", Options("variable=Precip parameters=testing/files/regression1order.txt"));
      EXPECT_TRUE(c);
      EXPECT_EQ("regression", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryQnh) {
      Calibrator* c;
      c = Calibrator::getScheme("qnh", Options("variable=Precip"));
      EXPECT_TRUE(c);
      EXPECT_EQ("qnh", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryValid) {
      Calibrator::getScheme("zaga", Options("variable=T parameters=testing/files/parameters.txt"));
      Calibrator::getScheme("zaga", Options("variable=Precip variable=T parameters=testing/files/parameters.txt"));
      Calibrator::getScheme("neighbourhood", Options("variable=Precip variable=T parameters=testing/files/parameters.txt"));
   }
   TEST_F(TestCalibrator, factoryInvalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(Calibrator::getScheme("zaga", Options("")), ".*");
      EXPECT_DEATH(Calibrator::getScheme("zaga", Options("parameters=testing/files/parameters.txt")), ".*");
      EXPECT_DEATH(Calibrator::getScheme("zaga", Options("variable=T")), ".*");
      EXPECT_DEATH(Calibrator::getScheme("cloud", Options("")), ".*");
      EXPECT_DEATH(Calibrator::getScheme("neighbourhood", Options("radius=-2")), ".*");
      EXPECT_DEATH(Calibrator::getScheme("phase", Options("")), ".*");
      EXPECT_DEATH(Calibrator::getScheme("regression", Options("")), ".*");
      EXPECT_DEATH(Calibrator::getScheme("regression", Options("variable=T")), ".*");
      EXPECT_DEATH(Calibrator::getScheme("regression", Options("parameters=testing/files/regression1order.txt")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
