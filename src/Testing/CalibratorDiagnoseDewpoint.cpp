#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/DiagnoseDewpoint.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorDiagnoseDewpoint : public ::testing::Test {
      protected:
         virtual void SetUp() {
            mT = Variable("air_temperature_2m");
            mRH = Variable("relativity_humidity_2m");
            mDewpoint = Variable("dewpoint");
            cal = new CalibratorDiagnoseDewpoint(mDewpoint ,Options("temperature=air_temperature_2m rh=relative_humidity_2m"));
         }
         virtual void TearDown() {
            delete cal;
         }
         Variable mT;
         Variable mRH;
         Variable mDewpoint;
         CalibratorDiagnoseDewpoint* cal;
   };
   TEST_F(TestCalibratorDiagnoseDewpoint, RH2dewpoint) {
      // EXPECT_FLOAT_EQ(273.15-1.27, CalibratorDiagnose::RH2dewpoint(273.15, 0.9));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnoseDewpoint::RH2dewpoint(Util::MV, 0.9));
      EXPECT_FLOAT_EQ(Util::MV, CalibratorDiagnoseDewpoint::RH2dewpoint(273.15, Util::MV));
   }
   TEST_F(TestCalibratorDiagnoseDewpoint, description) {
      CalibratorDiagnoseDewpoint::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
