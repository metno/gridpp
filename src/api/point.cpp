#include "gridpp.h"

using namespace gridpp;

gridpp::Point::Point(float lat, float lon, float elev, float laf, CoordinateType type) {
    this->lat = lat;
    this->lon = lon;
    this->elev = elev;
    this->laf = laf;
    this->type = type;

    if(type == gridpp::Geodetic) {
        float new_x = 0;
        float new_y = 0;
        float new_z = 0;
        gridpp::convert_coordinates(lat, lon, type, this->x, this->y, this->z);
    }
    else {
        this->x = lat;
        this->y = lon;
        this->z = 0;
    }
}

gridpp::Point::Point(float lat, float lon, float elev, float laf, CoordinateType type, float x, float y, float z) : lat(lat), lon(lon), elev(elev), laf(laf), type(type), x(x), y(y), z(z) {
}
