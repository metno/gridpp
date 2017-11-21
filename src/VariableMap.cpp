#include "VariableMap.h"
#include <sstream>
#include <fstream>
#include "Util.h"

VariableMap::VariableMap() {

}

void VariableMap::add(Variable::Type iType, std::string iName) {
   mMap[iType] = iName;
}

bool VariableMap::has(Variable::Type iType) const {
   std::map<int, std::string>::const_iterator it = mMap.find(iType);
   return it != mMap.end();
}

std::string VariableMap::get(Variable::Type iType) const {
   std::map<int, std::string>::const_iterator it = mMap.find(iType);
   if(it == mMap.end()) {
      std::stringstream ss;
      ss << "Cannot find variable " << Variable::getTypeName(iType) << " in map";
      Util::error(ss.str());
   }
   else {
      return it->second;
   }
}


bool VariableMap::load(std::string iFilename) {
   std::ifstream ifs(iFilename.c_str(), std::ifstream::in);
   if(!ifs.good()) {
      return false;
   }
   while(ifs.good()) {
      char line[10000];
      ifs.getline(line, 10000, '\n');
      if(ifs.good() && line[0] != '#') {
         std::stringstream ss(line);
         std::string key;
         std::string value;
         ss >> key;
         ss >> value;
         Variable::Type type = Variable::getType(key);
         add(type, value);
      }
   }
   return true;
}

void VariableMap::loadDefaults() {
   add(Variable::PrecipAcc, "precipitation_amount_acc");
   add(Variable::Cloud, "cloud_area_fraction");
   add(Variable::T, "air_temperature_2m");
   add(Variable::TMin, "air_temperature_2m_min6h");
   add(Variable::TMax, "air_temperature_2m_max6h");
   add(Variable::TD, "dew_point_temperature_2m");
   add(Variable::Tlevel0, "air_temperature_ml");
   add(Variable::Tlevel1, "air_temperature_ml");
   add(Variable::Precip, "precipitation_amount");
   add(Variable::Pop, "precipitation_amount_prob_low");
   add(Variable::Pop6h, "precipitation_amount_prob_low_6h");
   add(Variable::PrecipLow, "precipitation_amount_low_estimate");
   add(Variable::PrecipMiddle, "precipitation_amount_middle_estimate");
   add(Variable::PrecipHigh, "precipitation_amount_high_estimate");
   add(Variable::PrecipRate, "lwe_precipitation_rate");
   add(Variable::U, "eastward_wind_10m");
   add(Variable::Xwind, "x_wind_10m");
   add(Variable::V, "northward_wind_10m");
   add(Variable::Ywind, "y_wind_10m");
   add(Variable::W, "windspeed_10m");
   add(Variable::WD, "winddirection_10m");
   add(Variable::MSLP, "air_pressure_at_sea_level");
   add(Variable::RH, "relative_humidity_2m");
   add(Variable::Phase, "phase");
   add(Variable::P, "surface_air_pressure");
   add(Variable::MSLP, "air_pressure_at_sea_level");
   add(Variable::QNH, "qnh");
   add(Variable::SwinAcc, "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time");
   add(Variable::LwinAcc, "integral_of_surface_downwelling_longwave_flux_in_air_wrt_time");
   add(Variable::Fake, "fake");
}

void VariableMap::clear() {
   mMap.clear();
}
