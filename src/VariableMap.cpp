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


std::vector<Variable> VariableMap::load(std::string iFilename) {
   std::vector<Variable> variables;
   std::ifstream ifs(iFilename.c_str(), std::ifstream::in);
   if(!ifs.good()) {
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
            Variable variable(value, type);
            variables.push_back(variable);
         }
      }
   }
   return variables;
}

std::vector<Variable> VariableMap::loadDefaults() {
   std::vector<Variable> variables;
   variables.push_back(Variable("precipitation_amount_acc", Variable::PrecipAcc));
   variables.push_back(Variable("cloud_area_fraction", Variable::Cloud));
   variables.push_back(Variable("air_temperature_2m", Variable::T));
   variables.push_back(Variable("air_temperature_2m_min6h", Variable::TMin));
   variables.push_back(Variable("air_temperature_2m_max6h", Variable::TMax));
   variables.push_back(Variable("dew_point_temperature_2m", Variable::TD));
   variables.push_back(Variable("air_temperature_ml", Variable::Tlevel0));
   variables.push_back(Variable("air_temperature_ml", Variable::Tlevel1));
   variables.push_back(Variable("precipitation_amount", Variable::Precip));
   variables.push_back(Variable("precipitation_amount_prob_low", Variable::Pop));
   variables.push_back(Variable("precipitation_amount_prob_low_6h", Variable::Pop6h));
   variables.push_back(Variable("precipitation_amount_low_estimate", Variable::PrecipLow));
   variables.push_back(Variable("precipitation_amount_middle_estimate", Variable::PrecipMiddle));
   variables.push_back(Variable("precipitation_amount_high_estimate", Variable::PrecipHigh));
   variables.push_back(Variable("lwe_precipitation_rate", Variable::PrecipRate));
   variables.push_back(Variable("eastward_wind_10m", Variable::U));
   variables.push_back(Variable("x_wind_10m", Variable::Xwind));
   variables.push_back(Variable("northward_wind_10m", Variable::V));
   variables.push_back(Variable("y_wind_10m", Variable::Ywind));
   variables.push_back(Variable("windspeed_10m", Variable::W));
   variables.push_back(Variable("winddirection_10m", Variable::WD));
   variables.push_back(Variable("air_pressure_at_sea_level", Variable::MSLP));
   variables.push_back(Variable("relative_humidity_2m", Variable::RH));
   variables.push_back(Variable("phase", Variable::Phase));
   variables.push_back(Variable("surface_air_pressure", Variable::P));
   variables.push_back(Variable("air_pressure_at_sea_level", Variable::MSLP));
   variables.push_back(Variable("qnh", Variable::QNH));
   variables.push_back(Variable("integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time", Variable::SwinAcc));
   variables.push_back(Variable("integral_of_surface_downwelling_longwave_flux_in_air_wrt_time", Variable::LwinAcc));
   variables.push_back(Variable("fake", Variable::Fake));

   return variables;
}

static Variable getDefault(Variable::Type iType) {
   std::vector<Variable> variables = VariableMap::loadDefaults();
   for(int i = 0; i < variables.size(); i++) {
      if(iType == variables[i].getType())
         return variables[i];
   }
   // TODO:
   Util::error("Could not find variable of type");

}

void VariableMap::clear() {
   mMap.clear();
}
