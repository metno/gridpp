#include "gridpp.h" 
using namespace gridpp;
const float gridpp::StructureFunction::default_min_rho = 0.0013;
gridpp::StructureFunction::StructureFunction(float localization_distance) {
    if(!gridpp::is_valid(localization_distance) || localization_distance < 0)
        throw std::invalid_argument("Structure function initizlied with invalid localization distance");

    m_localization_distance = localization_distance;
}
float gridpp::StructureFunction::corr_background(const Point& p1, const Point& p2) const {
    return corr(p1, p2);
}
float gridpp::StructureFunction::barnes_rho(float dist, float length) const {
    if(!gridpp::is_valid(length) || length == 0)
        // Disabled
        return 1;
    if(!gridpp::is_valid(dist))
        return 0;
    float v = dist / length;
    return exp(-0.5 * v * v);
}
float gridpp::StructureFunction::cressman_rho(float dist, float length) const {
    if(!gridpp::is_valid(length) || length == 0)
        // Disabled
        return 1;
    if(!gridpp::is_valid(dist))
        return 0;
    if(dist >= length)
        return 0;
    return (length * length - dist * dist) / (length * length + dist * dist);
}
float gridpp::StructureFunction::localization_distance(const Point& p) const {
    return m_localization_distance;
}
gridpp::MultipleStructure::MultipleStructure(const StructureFunction& structure_h, const StructureFunction& structure_v, const StructureFunction& structure_w) {
    m_structure_h = structure_h.clone();
    m_structure_v = structure_v.clone();
    m_structure_w = structure_w.clone();
}
float gridpp::MultipleStructure::localization_distance(const Point& p) const {
    return m_structure_h->localization_distance(p);
}
float gridpp::MultipleStructure::corr(const Point& p1, const Point& p2) const {
    Point p1_h(p1.lat, p1.lon, p1.elev, p1.laf, p1.type);
    Point p2_h(p2.lat, p2.lon, p1.elev, p1.laf, p1.type);
    Point p1_v(p1.lat, p1.lon, p1.elev, p1.laf, p1.type);
    Point p2_v(p1.lat, p1.lon, p2.elev, p1.laf, p1.type);
    Point p1_w(p1.lat, p1.lon, p1.elev, p1.laf, p1.type);
    Point p2_w(p1.lat, p1.lon, p1.elev, p2.laf, p1.type);
    float corr_h = m_structure_h->corr(p1_h, p2_h);
    float corr_v = m_structure_v->corr(p1_v, p2_v);
    float corr_w = m_structure_w->corr(p1_w, p2_w);
    return corr_h * corr_v * corr_w;
}
gridpp::StructureFunction* gridpp::MultipleStructure::clone() const {
    gridpp::StructureFunction* val = new gridpp::MultipleStructure(*m_structure_h, *m_structure_v, *m_structure_w);
    return val;
}

/** Barnes */
gridpp::BarnesStructure::BarnesStructure(float h, float v, float w, float hmax) :
    gridpp::StructureFunction(h) {
    if(gridpp::is_valid(hmax) && hmax < 0)
        throw std::invalid_argument("hmax must be >= 0");
    if(!gridpp::is_valid(v) || v < 0)
        throw std::invalid_argument("v must be >= 0");
    if(!gridpp::is_valid(w) || w < 0)
        throw std::invalid_argument("w must be >= 0");

    // Override localization distance
    if(gridpp::is_valid(hmax))
        m_localization_distance = hmax;
    else {
        // Calculate the horizontal localization radius. For min_rho=0.0013, the factor is 3.64
        float default_min_rho = 0.0013;
        m_localization_distance = sqrt(-2*log(default_min_rho)) * h;
    }
    mH = h;
    mV = v;
    mW = w;
}
float gridpp::BarnesStructure::corr(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_distance_fast(p1, p2);
    if(hdist > localization_distance(p1))
        return 0;
    float rho = gridpp::StructureFunction::barnes_rho(hdist, mH);
    if(gridpp::is_valid(p1.elev) && gridpp::is_valid(p2.elev)) {
        float vdist = p1.elev - p2.elev;
        rho *= gridpp::StructureFunction::barnes_rho(vdist, mV);
    }
    if(gridpp::is_valid(p1.laf) && gridpp::is_valid(p2.laf)) {
        float lafdist = p1.laf - p2.laf;
        rho *= gridpp::StructureFunction::barnes_rho(lafdist, mW);
    }
    return rho;
}
gridpp::StructureFunction* gridpp::BarnesStructure::clone() const {
    gridpp::StructureFunction* val = new gridpp::BarnesStructure(mH, mV, mW, m_localization_distance);
    return val;
}

/** Cressman */
gridpp::CressmanStructure::CressmanStructure(float h, float v, float w) :
    gridpp::StructureFunction(h) {
    if(!gridpp::is_valid(v) || v < 0)
        throw std::invalid_argument("v must be >= 0");
    if(!gridpp::is_valid(w) || w < 0)
        throw std::invalid_argument("w must be >= 0");
    mH = h;
    mV = v;
    mW = w;
}
float gridpp::CressmanStructure::corr(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_distance_fast(p1, p2);
    float rho = gridpp::StructureFunction::cressman_rho(hdist, mH);
    if(gridpp::is_valid(p1.elev) && gridpp::is_valid(p2.elev)) {
        float vdist = p1.elev - p2.elev;
        rho *= gridpp::StructureFunction::cressman_rho(vdist, mV);
    }
    if(gridpp::is_valid(p1.laf) && gridpp::is_valid(p2.laf)) {
        float lafdist = p1.laf - p2.laf;
        rho *= gridpp::StructureFunction::cressman_rho(lafdist, mW);
    }
    return rho;
}
gridpp::StructureFunction* gridpp::CressmanStructure::clone() const {
    gridpp::StructureFunction* val = new gridpp::CressmanStructure(mH, mV, mW);
    return val;
}

/** Density */
gridpp::DensityStructure::DensityStructure(Grid grid, vec2 h, vec2 v, vec2 w, float min_rho) :
        m_grid(grid),
        m_min_rho(min_rho),
        mH(h),
        mV(v),
        mW(w) {
    if(grid.size()[0] != h.size() || grid.size()[0] != v.size() || grid.size()[0] != w.size())
        throw std::invalid_argument("Grid size not the same as scale size");
    if(grid.size()[1] != h[0].size() || grid.size()[1] != v[0].size() || grid.size()[1] != w[0].size())
        throw std::invalid_argument("Grid size not the same as scale size");
}
float gridpp::DensityStructure::corr(const Point& p1, const Point& p2) const {
    ivec I = m_grid.get_nearest_neighbour(p1.lat, p1.lon);
    if(I[0] > mH.size())
        throw std::runtime_error("Invalid I[0]");
    if(I[1] > mH[I[0]].size())
        throw std::runtime_error("Invalid I[1]");
    float h1 = mH[I[0]][I[1]];
    float v1 = mV[I[0]][I[1]];
    float w1 = mW[I[0]][I[1]];
#if 0
    I = m_grid.get_nearest_neighbour(p2.lat, p2.lon);
    if(I[0] > mH.size())
        throw std::runtime_error("Invalid i[0]");
    if(I[1] > mH[I[0]].size())
        throw std::runtime_error("Invalid I[1]");
    float h2 = mH[I[0]][I[1]];
    float v2 = mV[I[0]][I[1]];
    float w2 = mW[I[0]][I[1]];
    float h = (h1 + h2) / 2;
    float v = (v1 + v2) / 2;
    float w = (w1 + w2) / 2;
#else
    float h = h1;
    float v = v1;
    float w = w1;
#endif

    float hdist = gridpp::KDTree::calc_distance_fast(p1, p2);
    // TODO:
    if(hdist > localization_distance(p1))
        return 0;
    float rho = gridpp::StructureFunction::barnes_rho(hdist, h);
    if(gridpp::is_valid(p1.elev) && gridpp::is_valid(p2.elev)) {
        float vdist = p1.elev - p2.elev;
        rho *= gridpp::StructureFunction::barnes_rho(vdist, v);
    }
    if(gridpp::is_valid(p1.laf) && gridpp::is_valid(p2.laf)) {
        float lafdist = p1.laf - p2.laf;
        rho *= gridpp::StructureFunction::barnes_rho(lafdist, w);
    }
    return rho;
}
float gridpp::DensityStructure::corr_background(const Point& p1, const Point& p2) const {
    return corr(p1, p2);
    float hdist = gridpp::KDTree::calc_distance_fast(p1, p2);
    // TODO:
    if(hdist > localization_distance(p1))
        return 0;
    // TODO
    float rho = gridpp::StructureFunction::barnes_rho(hdist, 10000);
    return rho;
}
float gridpp::DensityStructure::localization_distance(const Point& p) const {
    ivec I = m_grid.get_nearest_neighbour(p.lat, p.lon);
    float curr = sqrt(-2*log(m_min_rho)) * mH[I[0]][I[1]];
    return curr;
}
gridpp::StructureFunction* gridpp::DensityStructure::clone() const {
    gridpp::StructureFunction* val = new gridpp::DensityStructure(m_grid, mH, mV, mW, m_min_rho);
    return val;
}

/** CrossValidation */
gridpp::CrossValidation::CrossValidation(StructureFunction& structure, float dist) :
        StructureFunction(0){
    if(!gridpp::is_valid(dist) || dist < 0)
        throw std::invalid_argument("Invalid 'dist' in CrossValidation structure");
    m_structure = structure.clone();
    m_dist = dist;
}
float gridpp::CrossValidation::corr(const Point& p1, const Point& p2) const {
    return m_structure->corr(p1, p2);
}
float gridpp::CrossValidation::corr_background(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_distance_fast(p1, p2);
    if(gridpp::is_valid(m_dist)) {
        if(hdist <= m_dist)
            return 0;
    }
    return m_structure->corr_background(p1, p2);
}
gridpp::StructureFunction* gridpp::CrossValidation::clone() const {
    gridpp::StructureFunction* val = new gridpp::CrossValidation(*m_structure);
    return val;
}
float gridpp::CrossValidation::localization_distance(const Point& p) const {
    return m_structure->localization_distance(p);
}
