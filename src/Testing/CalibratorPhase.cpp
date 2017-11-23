#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Phase.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorPhase : public ::testing::Test {
      protected:
         TestCalibratorPhase() {
         }
         virtual ~TestCalibratorPhase() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
         Parameters getParameters(float a1, float a2) {
            std::vector<float> parValues(8, 0);
            parValues[0] = a1;
            parValues[1] = a2;
            return Parameters (parValues);
         }
         ParameterFileSimple getParameterFile(float a1, float a2) {
            std::vector<float> values;
            values.push_back(a1);
            values.push_back(a2);
            ParameterFileSimple parFile(values);
            return parFile;
         }
         CalibratorPhase getCalibratorWetbulb(Variable iVariable) {
            return CalibratorPhase(iVariable, Options("temperatureVariable=air_temperature_2m precipitationVariable=precipitation_amount rhVariable=relative_humidity_2m pressureVariable=surface_air_pressure"));
         }
         CalibratorPhase getCalibrator(Variable iVariable) {
            return CalibratorPhase(iVariable, Options("temperatureVariable=air_temperature_2m precipitationVariable=precipitation_amount"));
         }
         void setValues(const File& iFile, float iPrecip, float iTemp, float iRh, float iPressure) {
            FieldPtr precip   = iFile.getField(Variable("precipitation_amount"), 0);
            FieldPtr temp     = iFile.getField(Variable("air_temperature_2m"), 0);
            FieldPtr rh       = iFile.getField(Variable("relative_humidity_2m"), 0);
            FieldPtr pressure = iFile.getField(Variable("surface_air_pressure"), 0);
            (*precip)(0,0,0) = iPrecip;
            (*temp)(0,0,0) = iTemp;
            (*rh)(0,0,0) = iRh;
            (*pressure)(0,0,0) = iPressure;
         }
   };

   TEST_F(TestCalibratorPhase, 10x10) {
      FileNetcdf file("testing/files/10x10.nc");
      ParameterFileText parFile(Options("file=testing/files/parametersPhase.txt"));
      Parameters parameters = parFile.getParameters(0);
      ASSERT_EQ(2, parameters.size());
      EXPECT_FLOAT_EQ(273.7, parameters[0]);
      EXPECT_FLOAT_EQ(274.7, parameters[1]);
      CalibratorPhase cal = getCalibratorWetbulb(Variable("phase"));

      cal.calibrate(file, &parFile);
      FieldPtr phase    = file.getField(Variable("phase"), 0);
      FieldPtr precip   = file.getField(Variable("precipitation_amount"), 0);
      FieldPtr temp     = file.getField(Variable("air_temperature_2m"), 0);
      FieldPtr rh       = file.getField(Variable("relative_humidity_2m"), 0);
      FieldPtr pressure = file.getField(Variable("surface_air_pressure"), 0);
      // T      301 K
      // RH     0.95925
      // P      98334 pa
      // Precip 1.38 mm
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseRain, (*phase)(2,5,0));
   }
   TEST_F(TestCalibratorPhase, phases) {
      FileFake file(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      //file.initNewVariable(Variable("phase"));
      //file.initNewVariable(Variable("precipitation_amount"));
      //file.initNewVariable(Variable("air_temperature_2m"));
      //file.initNewVariable(Variable("relative_humidity_2m"));
      //file.initNewVariable(Variable("surface_air_pressure"));
      FieldPtr phase = file.getField(Variable("phase"), 0);
      ParameterFileSimple parFile = getParameterFile(273.7,274.7);
      CalibratorPhase cal = getCalibratorWetbulb(Variable("phase"));
      setValues(file, 5, 270, 0.95, 101325);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSnow, (*phase)(0,0,0));

      setValues(file, 5, 274, 0.98, 101325);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSleet, (*phase)(0,0,0));

      setValues(file, 0, 274, 0.5, 101325);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseNone, (*phase)(0,0,0));

      setValues(file, 5, 275, 0.5, 101325);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSnow, (*phase)(0,0,0));
   }
   TEST_F(TestCalibratorPhase, useWetbulb) {
      FileFake file(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      FieldPtr phase = file.getField(Variable("phase"), 0);
      ParameterFileSimple parFile = getParameterFile(273.7,274.7);
      CalibratorPhase calwet = getCalibratorWetbulb(Variable("phase"));
      CalibratorPhase cal = getCalibrator(Variable("phase"));
      setValues(file, 5, 275, 0.5, 101325);

      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseRain, (*phase)(0,0,0));

      calwet.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSnow, (*phase)(0,0,0));

      setValues(file, 5, 274, 0.5, 101325);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSleet, (*phase)(0,0,0));
      calwet.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSnow, (*phase)(0,0,0));
   }
   TEST_F(TestCalibratorPhase, missingParameters) {
      FileNetcdf file("testing/files/10x10.nc");
      ParameterFileSimple parFile = getParameterFile(Util::MV,1.5);
      CalibratorPhase cal = getCalibrator(Variable("phase"));

      cal.calibrate(file, &parFile);
      FieldPtr phase = file.getField(Variable("phase"), 0);
      for(int i = 0; i < file.getNumLat(); i++) {
         for(int j = 0; j < file.getNumLon(); j++) {
            for(int e = 0; e < file.getNumEns(); e++) {
               EXPECT_FLOAT_EQ(Util::MV, (*phase)(i,j,e));
            }
         }
      }
   }
   TEST_F(TestCalibratorPhase, getWetbulb) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
                                                                                  // NOAA
      EXPECT_FLOAT_EQ(269.02487, CalibratorPhase::getWetbulb(270, 100000, 0.80)); // 269.03
      EXPECT_FLOAT_EQ(296.13763, CalibratorPhase::getWetbulb(300, 101000, 0.70)); // 295.95
      EXPECT_FLOAT_EQ(269.92218, CalibratorPhase::getWetbulb(270, 100000, 1));    // 270
      EXPECT_FLOAT_EQ(239.83798, CalibratorPhase::getWetbulb(240, 50000, 0.90));  // 239.89
   }
   TEST_F(TestCalibratorPhase, getWetbulbInvalid) {
      EXPECT_FLOAT_EQ(Util::MV, CalibratorPhase::getWetbulb(Util::MV, 30000, 1));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorPhase::getWetbulb(270, Util::MV, 1));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorPhase::getWetbulb(270, 30000, Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorPhase::getWetbulb(270, 100000, 0)); // No humidity
   }
   TEST_F(TestCalibratorPhase, description) {
      CalibratorPhase::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
