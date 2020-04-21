#include "gridpp.h"

gridpp::Parameters::Parameters(const gridpp::Points& points, const vec2& values, gridpp::Interpolator& interpolator) :
    m_points(points),
    m_values(values)
    {
    m_interpolator = interpolator.clone();
}
gridpp::Parameters::~Parameters() {
    delete m_interpolator;
}

vec gridpp::Parameters::get(float lat, float lon, float altitude, float land_area_fraction) const {
    int P = m_values.size();
    vec values(P);
    for(int p = 0; p < P; p++) {
        values[p] = m_interpolator->interpolate(m_points, m_values[p], lat, lon, altitude, land_area_fraction);
    }
    return values;
}
