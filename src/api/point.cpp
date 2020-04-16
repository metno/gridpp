#include "gridpp.h"

gridpp::Point::Point(float lat, float lon, float elev, float laf) {
    this->lat = lat;
    this->lon = lon;
    this->elev = elev;
    this->laf = laf;
}
