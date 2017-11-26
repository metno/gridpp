#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/DiagnoseHumidity.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorDiagnoseHumidity : public ::testing::Test {
      protected:
         virtual void SetUp() {
            mTemperature = Variable("air_temperature_2m");
            mRh = Variable("relative_humidity_2m");
            mDewpoint = Variable("dewpoint");
            mWetbulb = Variable("wetbulb");
            mPressure = Variable("surface_air_pressure");
            calRh = new CalibratorDiagnoseHumidity(mRh, Options("temperature=air_temperature_2m dewpoint=dewpoint compute=rh"));
            calDewpoint = new CalibratorDiagnoseHumidity(mDewpoint, Options("temperature=air_temperature_2m rh=relative_humidity_2m compute=dewpoint"));
            calWetbulb  = new CalibratorDiagnoseHumidity(mWetbulb, Options("temperature=air_temperature_2m rh=relative_humidity_2m compute=wetbulb"));
            calWetbulbP = new CalibratorDiagnoseHumidity(mWetbulb, Options("temperature=air_temperature_2m rh=relative_humidity_2m pressure=surface_air_pressure compute=wetbulb"));
         }
         virtual void TearDown() {
            delete calRh;
            delete calDewpoint;
            delete calWetbulb;
            delete calWetbulbP;
         }
         Variable mTemperature;
         Variable mRh;
         Variable mDewpoint;
         Variable mPressure;
         Variable mWetbulb;
         CalibratorDiagnoseHumidity* calRh;
         CalibratorDiagnoseHumidity* calDewpoint;
         CalibratorDiagnoseHumidity* calWetbulb;
         CalibratorDiagnoseHumidity* calWetbulbP;
   };
   TEST_F(TestCalibratorDiagnoseHumidity, dewpoint) {
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnoseHumidity::computeDewpoint(Util::MV, 0.9));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnoseHumidity::computeDewpoint(273.15, Util::MV));
   }
   TEST_F(TestCalibratorDiagnoseHumidity, rh) {
      // EXPECT_FLOAT_EQ(0.9, CalibratorDiagnose::dewpoint2RH(273.15, 273.15-1.27));
      EXPECT_FLOAT_EQ(1, CalibratorDiagnoseHumidity::computeRh(273.15, 273.15));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnoseHumidity::computeRh(Util::MV, 273.15-1.27));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnoseHumidity::computeRh(273.15, Util::MV));
   }
   TEST_F(TestCalibratorDiagnoseHumidity, description) {
      CalibratorDiagnoseHumidity::description();
   }
   TEST_F(TestCalibratorDiagnoseHumidity, valid) {
      CalibratorDiagnoseHumidity(mRh, Options("temperature=air_temperature_2m dewpoint=dewpoint compute=rh"));
      CalibratorDiagnoseHumidity(mDewpoint, Options("temperature=air_temperature_2m rh=relative_humidity_2m compute=dewpoint"));
      CalibratorDiagnoseHumidity(mWetbulb, Options("temperature=air_temperature_2m rh=relative_humidity_2m compute=wetbulb"));
      CalibratorDiagnoseHumidity(mWetbulb, Options("temperature=air_temperature_2m rh=relative_humidity_2m pressure=surface_air_pressure compute=wetbulb"));
   }
   TEST_F(TestCalibratorDiagnoseHumidity, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // Missing variables for RH
      EXPECT_DEATH(CalibratorDiagnoseHumidity(mRh, Options("temperature=air_temperature_2m compute=rh")), ".*");
      EXPECT_DEATH(CalibratorDiagnoseHumidity(mRh, Options("temperature=air_temperature_2m pressure=pressure rh=relative_humidity compute=rh")), ".*");
      // Missing variables for dewpoint
      EXPECT_DEATH(CalibratorDiagnoseHumidity(mDewpoint, Options("temperature=air_temperature_2m compute=dewpoint")), ".*");
      EXPECT_DEATH(CalibratorDiagnoseHumidity(mDewpoint, Options("temperature=air_temperature_2m pressure=pressure dewpoint=dewpoint compute=dewpoint")), ".*");
      // Missing variables for wetbulb
      EXPECT_DEATH(CalibratorDiagnoseHumidity(mDewpoint, Options("temperature=air_temperature_2m compute=wetbulb")), ".*");
      EXPECT_DEATH(CalibratorDiagnoseHumidity(mDewpoint, Options("temperature=air_temperature_2m pressure=pressure dewpoint=dewpoint compute=wetbulb")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
