#include "gridpp.h"

gridpp::StructureFunction::StructureFunction(float localization_distance) {
    if(!gridpp::util::is_valid(localization_distance) || localization_distance < 0)
        throw std::invalid_argument("Invalid 'localization_distance' in structure");

    mLocalizationDistance = localization_distance;
}
float gridpp::StructureFunction::corr(const Point& p1, const Point& p2, const vec& v1, const vec& v2) const {
    return corr(p1, p2);
}
float gridpp::StructureFunction::barnes_rho(float dist, float length) const {
    float v = dist / length;
    return exp(-0.5 * v * v);
}
float gridpp::StructureFunction::barnes_rho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wlength) const {
    float h = (iHDist/hlength);
    float rho = exp(-0.5 * h * h);
    if(gridpp::util::is_valid(vlength) && vlength > 0) {
        if(!gridpp::util::is_valid(iVDist)) {
            rho = 0;
        }
        else {
            float v = (iVDist/vlength);
            float factor = exp(-0.5 * v * v);
            rho *= factor;
        }
    }
    if(gridpp::util::is_valid(wlength) && wlength > 0) {
        float factor = exp(-0.5 * iLDist * iLDist / (wlength * wlength));
        rho *= factor;
    }
    return rho;
}
float gridpp::StructureFunction::cressman_rho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wlength) const {
    float h = (iHDist/hlength);
    float rho = (hlength*hlength - iHDist * iHDist) / (hlength*hlength + iHDist * iHDist);
    if(iHDist > hlength)
        rho = 0;
    if(gridpp::util::is_valid(vlength) && vlength > 0) {
        if(!gridpp::util::is_valid(iVDist)) {
            rho = 0;
        }
        else {
            float factor = (vlength*vlength - iVDist * iVDist) / (vlength*vlength + iVDist * iVDist);
            rho *= factor;
        }
    }
    if(gridpp::util::is_valid(wlength) && wlength > 0) {
        float factor = (wlength*wlength - iLDist * iLDist) / (wlength*wlength + iLDist * iLDist);
        rho *= factor;
    }
    return rho;
}
float gridpp::StructureFunction::cressman_rho(float dist, float length) const {
    return (length * length - dist * dist) / (length * length + dist * dist);
}
float gridpp::StructureFunction::localization_distance() const {
    return mLocalizationDistance;
}
gridpp::BarnesStructure::BarnesStructure(float h, float v, float w, float hmax) :
    gridpp::StructureFunction(h) {
    if(gridpp::util::is_valid(hmax))
        mLocalizationDistance = hmax;
    else {
        // Calculate the horizontal localization radius. For min_rho=0.0013, the factor is 3.64
        float default_min_rho = 0.0013;
        mLocalizationDistance = sqrt(-2*log(default_min_rho)) * h;
    }
    mH = h;
    mV = v;
    mW = w;
}
float gridpp::BarnesStructure::corr(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_distance(p1.lat, p1.lon, p2.lat, p2.lon);
    if(hdist > localization_distance())
        return 0;
    float vdist = gridpp::MV;
    if(gridpp::util::is_valid(p1.elev) && gridpp::util::is_valid(p2.elev))
        vdist = p1.elev - p2.elev;
    float lafdist = 0;
    if(gridpp::util::is_valid(p1.laf) && gridpp::util::is_valid(p2.laf))
        lafdist = p1.laf - p2.laf;
    float rho = barnes_rho(hdist, vdist, lafdist, mH, mV, mW);
    return rho;
}
gridpp::CressmanStructure::CressmanStructure(float h, float v, float w) :
    gridpp::StructureFunction(h) {
    mH = h;
    mV = v;
    mW = w;
}
float gridpp::CressmanStructure::corr(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_distance(p1.lat, p1.lon, p2.lat, p2.lon);
    float vdist = gridpp::MV;
    if(gridpp::util::is_valid(p1.elev) && gridpp::util::is_valid(p2.elev))
        vdist = p1.elev - p2.elev;
    float lafdist = 0;
    if(gridpp::util::is_valid(p1.laf) && gridpp::util::is_valid(p2.laf))
        lafdist = p1.laf - p2.laf;
    float rho = cressman_rho(hdist, vdist, lafdist, mH, mV, mW);
    return rho;
}
