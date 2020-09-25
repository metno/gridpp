#include "gridpp.h"

using namespace gridpp;

gridpp::Point::Point(float lat, float lon, float elev, float laf, bool flat) {
    this->lat = lat;
    this->lon = lon;
    this->elev = elev;
    this->laf = laf;
    this->flat = flat;
}
