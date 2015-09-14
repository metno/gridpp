#include "../ParameterFile/ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   TEST(ParameterFileNetcdfTest, singleTime) {
      ParameterFileNetcdf file(Options("file=testing/files/10x10_param.nc"));
      ASSERT_EQ(2, file.getTimes().size());
      ASSERT_EQ(100, file.getLocations().size());
      // Location loc(5,5,142.1456);
      Location loc(5,5,3);
      Parameters par = file.getParameters(0, loc);
      ASSERT_EQ(2, par.size());
      EXPECT_FLOAT_EQ(0.3, par[0]);
      EXPECT_FLOAT_EQ(2.3, par[1]);
   }
   TEST(ParameterFileNetcdfTest, singleTime_xy) {
      ParameterFileNetcdf file(Options("file=testing/files/10x10_param_xy.nc"));
      ASSERT_EQ(2, file.getTimes().size());
      ASSERT_EQ(100, file.getLocations().size());
      // Location loc(5,5,142.1456);
      Location loc(5,5,3);
      Parameters par = file.getParameters(0, loc);
      ASSERT_EQ(2, par.size());
      EXPECT_FLOAT_EQ(0.3, par[0]);
      EXPECT_FLOAT_EQ(2.3, par[1]);
   }
   /*
   // 10x10_param_xy.nc has longitude and coefficients with x,y ordering
   // Should still give you the same results
   TEST(ParameterFileNetcdfTest, xy_order) {
      ParameterFileNetcdf file_yx(Options("file=testing/files/10x10_param.nc"));
      ParameterFileNetcdf file_xy(Options("file=testing/files/10x10_param_xy.nc"));
      EXPECT_EQ(file_yx.getTimes(), file_xy.getTimes());
      std::vector<Location> loc_yx = file_yx.getLocations();
      std::vector<Location> loc_xy = file_yx.getLocations();
      for(int i = 0; i < loc_yx.size(); i++) {
         if(loc_yx[i].lat() == 5 && loc_yx[i].lon() == 5) {
            Parameters par_yx= file_yx.getParameters(0, loc_yx[i]);
            Parameters par_xy= file_xy.getParameters(0, loc_yx[i]);
            std::cout << i << " " << loc_yx[i].lat() << " " << loc_yx[i].lon() << std::endl;
            std::cout << i << " " << loc_xy[i].lat() << " " << loc_xy[i].lon() << std::endl;
            EXPECT_EQ(par_yx.getValues(), par_xy.getValues());
            std::cout << par_yx[0] << std::endl;
            std::cout << par_xy[0] << std::endl;
         }
      }
   }
   */
   TEST(ParameterFileNetcdfTest, description) {
      ParameterFileNetcdf::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
