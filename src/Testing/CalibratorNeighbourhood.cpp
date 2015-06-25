#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Neighbourhood.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorNeighbourhood : public ::testing::Test {
      protected:
         TestCalibratorNeighbourhood() {
         }
         virtual ~TestCalibratorNeighbourhood() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
   };
   TEST_F(TestCalibratorNeighbourhood, 10x10_double) {
      FileArome from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T ,Options("radius=1"));

      cal.calibrate(from);
      FieldPtr after = from.getField(Variable::T, 0);
      ASSERT_EQ(10, after->getNumLat());
      ASSERT_EQ(10, after->getNumLon());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(304.6667, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(306.1667, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(308.25,   (*after)(0,9,0));

      CalibratorNeighbourhood cal2 = CalibratorNeighbourhood(Variable::T ,Options("radius=2"));
      cal2.calibrate(from);
      EXPECT_FLOAT_EQ(304.73114, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(305.35556, (*after)(5,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, 10x10_missingValues) {
      FileArome from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T, Options("radius=1"));
      FieldPtr field = from.getField(Variable::T, 0);
      (*field)(4,1,0) = Util::MV;
      (*field)(4,2,0) = Util::MV;
      (*field)(4,3,0) = Util::MV;
      (*field)(5,1,0) = Util::MV;
      (*field)(5,2,0) = Util::MV;
      (*field)(5,3,0) = Util::MV;
      (*field)(6,1,0) = Util::MV;
      (*field)(6,2,0) = Util::MV;
      (*field)(6,3,0) = Util::MV;

      cal.calibrate(from);
      EXPECT_FLOAT_EQ(Util::MV, (*field)(5,2,0));
      EXPECT_FLOAT_EQ(304.6667, (*field)(5,1,0));
   }
   TEST_F(TestCalibratorNeighbourhood, getRadius) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T, Options(""));
      // Check that default is valid
      EXPECT_GE(cal.getRadius(), 1);

      CalibratorNeighbourhood cal12 = CalibratorNeighbourhood(Variable::T, Options("radius=12"));
      EXPECT_FLOAT_EQ(12, cal12.getRadius());
   }
   TEST_F(TestCalibratorNeighbourhood, mean) {
      FileArome from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T ,Options("radius=1 stat=mean"));

      cal.calibrate(from);
      FieldPtr after = from.getField(Variable::T, 0);

      EXPECT_FLOAT_EQ(304.6667, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(306.1667, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(308.25,   (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, min) {
      FileArome from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T ,Options("radius=1 stat=min"));

      cal.calibrate(from);
      FieldPtr after = from.getField(Variable::T, 0);

      EXPECT_FLOAT_EQ(301, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(302, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(302, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, max) {
      FileArome from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T ,Options("radius=1 stat=max"));

      cal.calibrate(from);
      FieldPtr after = from.getField(Variable::T, 0);

      EXPECT_FLOAT_EQ(310, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(310, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(320, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, std) {
      FileArome from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T ,Options("radius=1 stat=std"));

      cal.calibrate(from);
      FieldPtr after = from.getField(Variable::T, 0);

      EXPECT_FLOAT_EQ(2.7080128, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(3.2871804, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(7.0133801, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, median) {
      FileArome from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T ,Options("radius=1 stat=median"));

      cal.calibrate(from);
      FieldPtr after = from.getField(Variable::T, 0);

      EXPECT_FLOAT_EQ(306, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(306, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(305.5, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, quantile) {
      FileArome from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(Variable::T ,Options("radius=1 stat=quantile quantile=0.9"));

      cal.calibrate(from);
      FieldPtr after = from.getField(Variable::T, 0);

      EXPECT_FLOAT_EQ(306.8, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(310,   (*after)(5,9,0));
      EXPECT_FLOAT_EQ(316.1, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      EXPECT_DEATH(CalibratorNeighbourhood(Variable::T, Options("stat=invalidStatType")), ".*");

      // Negative radius
      EXPECT_DEATH(CalibratorNeighbourhood(Variable::T, Options("radius=-1")), ".*");
      EXPECT_DEATH(CalibratorNeighbourhood(Variable::T, Options("radius=-1 stat=quantile")), ".*");
      EXPECT_DEATH(CalibratorNeighbourhood(Variable::T, Options("radius=-1 stat=quantile quantile=0.5")), ".*");

      // Missing quantile
      EXPECT_DEATH(CalibratorNeighbourhood(Variable::T, Options("stat=quantile")), ".*");
      // Invalid quantile
      EXPECT_DEATH(CalibratorNeighbourhood(Variable::T, Options("stat=quantile quantile=-0.1")), ".*");
      EXPECT_DEATH(CalibratorNeighbourhood(Variable::T, Options("stat=quantile quantile=1.1")), ".*");
      EXPECT_DEATH(CalibratorNeighbourhood(Variable::T, Options("stat=quantile quantile=-999")), ".*");
   }
   TEST_F(TestCalibratorNeighbourhood, description) {
      CalibratorNeighbourhood::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
