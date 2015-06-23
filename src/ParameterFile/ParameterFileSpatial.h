#ifndef PARAMETER_FILE_SPATIAL_H
#define PARAMETER_FILE_SPATIAL_H
#include <iostream>
#include <map>
#include "ParameterFile.h"
#include "../Parameters.h"
#include "../Location.h"

//! Represents parameters for different spatial coordinates.
class ParameterFileSpatial {
   public:
      virtual Parameters getParameters(int iTime, const Location& iLocation) const = 0;
      virtual std::vector<Location> getLocations() const = 0;
      virtual std::vector<int> getTimes() const = 0;
};
#include "ParameterFileSpatialSimple.h"
#include "ParameterFileSpatialMetnoKalman.h"
#endif
