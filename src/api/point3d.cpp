#include "gridpp.h"

using namespace gridpp;

gridpp::Point3D::Point3D(float x, float y, float z, float lat, float lon, float elev, float laf, CoordinateType type) : x(x), y(y), z(z), lat(lat), lon(lon), elev(elev), laf(laf), type(type)  {
}

gridpp::Point3D::Point3D(float lat, float lon, float elev, float laf, CoordinateType type) : x(0), y(0), z(0), lat(lat), lon(lon), elev(elev), laf(laf), type(type)  {
    if(type == gridpp::Geodetic) {
        float new_x = 0;
        float new_y = 0;
        float new_z = 0;
        gridpp::convert_coordinates(lat, lon, new_x, new_y, new_z, type);
        set(new_x, new_y, new_z, lat, lon, elev, laf, type);
    }
    else {
        set(lat, lon, 0, lat, lon, elev, laf, type);
    }
}
void gridpp::Point3D::set(float x, float y, float z, float lat, float lon, float elev, float laf, CoordinateType type) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->lat = lat;
    this->lon = lon;
    this->elev = elev;
    this->laf = laf;
    this->type = type;
}
