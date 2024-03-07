#include "gridpp.h"

using namespace gridpp;

const float gridpp::StructureFunction::default_min_rho = 0.0013;

gridpp::StructureFunction::StructureFunction(float localization_distance) {
    if(!gridpp::is_valid(localization_distance) || localization_distance < 0)
        throw std::invalid_argument("Structure function initizlied with invalid localization distance");

    m_localization_distance = localization_distance;
}
vec gridpp::StructureFunction::corr(const Point& p1, const std::vector<Point>& p2) const {
    vec output(p2.size());
    for(int i = 0; i < p2.size(); i++) {
        output[i] = corr(p1, p2[i]);
    }
    return output;
}
float gridpp::StructureFunction::corr_background(const Point& p1, const Point& p2) const {
    return corr(p1, p2);
}
vec gridpp::StructureFunction::corr_background(const Point& p1, const std::vector<Point>& p2) const {
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
    // Perturn coordinate only
    Point p1_h(p1.lat, p1.lon, p1.elev, p1.laf, p1.type, p1.x, p1.y, p1.z);
    Point p2_h(p2.lat, p2.lon, p1.elev, p1.laf, p1.type, p2.x, p2.y, p2.z);
    // Perturb elev only
    Point p1_v(p1.lat, p1.lon, p1.elev, p1.laf, p1.type, p1.x, p1.y, p1.z);
    Point p2_v(p1.lat, p1.lon, p2.elev, p1.laf, p1.type, p1.x, p1.y, p1.z);
    // Perturb laf only
    Point p1_w(p1.lat, p1.lon, p1.elev, p1.laf, p1.type, p1.x, p1.y, p1.z);
    Point p2_w(p1.lat, p1.lon, p1.elev, p2.laf, p1.type, p1.x, p1.y, p1.z);
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
    m_is_spatial(false) {
    if(gridpp::is_valid(hmax) && hmax < 0)
        throw std::invalid_argument("hmax must be >= 0");
    if(!gridpp::is_valid(h) || h < 0)
        throw std::invalid_argument("h must be >= 0");
    if(!gridpp::is_valid(v) || v < 0)
        throw std::invalid_argument("v must be >= 0");
    if(!gridpp::is_valid(w) || w < 0)
        throw std::invalid_argument("w must be >= 0");

    if(gridpp::is_valid(hmax))
        m_min_rho = exp(pow(hmax / h, 2)/-2);
    else
        m_min_rho = default_min_rho;
    vec2 h2(1);
    h2[0].push_back(h);
    vec2 v2(1);
    v2[0].push_back(v);
    vec2 w2(1);
    w2[0].push_back(w);
    mH = h2;
    mV = v2;
    mW = w2;
}
gridpp::BarnesStructure::BarnesStructure(Grid grid, vec2 h, vec2 v, vec2 w, float min_rho) :
        m_grid(grid),
        m_min_rho(min_rho),
        mH(h),
        mV(v),
        mW(w) {
    if(mH.size() == 1 && mH[0].size() == 1 && mV.size() == 1 && mV[0].size() == 1 && mW.size() == 1 && mW[0].size() == 1) {
        m_is_spatial = false;
    }
    else {
        m_is_spatial = true;
        if(grid.size()[0] != h.size() || grid.size()[0] != v.size() || grid.size()[0] != w.size())
            throw std::invalid_argument("Grid size not the same as scale size");
        if(grid.size()[1] != h[0].size() || grid.size()[1] != v[0].size() || grid.size()[1] != w[0].size())
            throw std::invalid_argument("Grid size not the same as scale size");
    }
}
float gridpp::BarnesStructure::corr(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_straight_distance(p1, p2);
    float rho = 1;
    if(m_is_spatial) {
        // This part is slower because of the nearest neighbour lookup
        ivec I = m_grid.get_nearest_neighbour(p1.lat, p1.lon);
        if(I[0] > mH.size())
            throw std::runtime_error("Invalid I[0]");
        if(I[1] > mH[I[0]].size())
            throw std::runtime_error("Invalid I[1]");

        float h = mH[I[0]][I[1]];
        float v = mV[I[0]][I[1]];
        float w = mW[I[0]][I[1]];

        // Use h to compute this, so we don't have to get nearest neighbour twice
        float loc_dist = localization_distance(h);
        if(hdist > loc_dist)
            return 0;

        rho = gridpp::StructureFunction::barnes_rho(hdist, h);
        if(gridpp::is_valid(p1.elev) && gridpp::is_valid(p2.elev)) {
            float vdist = p1.elev - p2.elev;
            rho *= gridpp::StructureFunction::barnes_rho(vdist, v);
        }
        if(gridpp::is_valid(p1.laf) && gridpp::is_valid(p2.laf)) {
            float lafdist = p1.laf - p2.laf;
            rho *= gridpp::StructureFunction::barnes_rho(lafdist, w);
        }
    }
    else {
        if(hdist > localization_distance(p1))
            return 0;

        rho = gridpp::StructureFunction::barnes_rho(hdist, mH[0][0]);
        if(gridpp::is_valid(p1.elev) && gridpp::is_valid(p2.elev)) {
            float vdist = p1.elev - p2.elev;
            rho *= gridpp::StructureFunction::barnes_rho(vdist, mV[0][0]);
        }
        if(gridpp::is_valid(p1.laf) && gridpp::is_valid(p2.laf)) {
            float lafdist = p1.laf - p2.laf;
            rho *= gridpp::StructureFunction::barnes_rho(lafdist, mW[0][0]);
        }
    }
    return rho;
}
vec gridpp::BarnesStructure::corr(const Point& p1, const std::vector<Point>& p2) const {
    vec output(p2.size());
    if(m_is_spatial) {
        ivec I = m_grid.get_nearest_neighbour(p1.lat, p1.lon);
        if(I[0] > mH.size())
            throw std::runtime_error("Invalid I[0]");
        if(I[1] > mH[I[0]].size())
            throw std::runtime_error("Invalid I[1]");

        float h = mH[I[0]][I[1]];
        float v = mV[I[0]][I[1]];
        float w = mW[I[0]][I[1]];
        float loc_dist = localization_distance(h);
        for(int i = 0; i < p2.size(); i++) {
            float hdist = gridpp::KDTree::calc_straight_distance(p1, p2[i]);
            float rho = 0;
            if(hdist <= loc_dist) {
                rho = gridpp::StructureFunction::barnes_rho(hdist, h);
                if(gridpp::is_valid(p1.elev) && gridpp::is_valid(p2[i].elev)) {
                    float vdist = p1.elev - p2[i].elev;
                    rho *= gridpp::StructureFunction::barnes_rho(vdist, v);
                }
                if(gridpp::is_valid(p1.laf) && gridpp::is_valid(p2[i].laf)) {
                    float lafdist = p1.laf - p2[i].laf;
                    rho *= gridpp::StructureFunction::barnes_rho(lafdist, w);
                }
            }
            output[i] = rho;
        }
    }
    else {
        for(int i = 0; i < p2.size(); i++) {
            output[i] = corr(p1, p2[i]);
        }
    }
    return output;
}
gridpp::StructureFunction* gridpp::BarnesStructure::clone() const {
    gridpp::StructureFunction* val = new gridpp::BarnesStructure(m_grid, mH, mV, mW, m_min_rho);
    return val;
}
float gridpp::BarnesStructure::localization_distance(const Point& p) const {
    if(m_is_spatial) {
        ivec I = m_grid.get_nearest_neighbour(p.lat, p.lon);
        return localization_distance(mH[I[0]][I[1]]);
    }
    else {
        return localization_distance(mH[0][0]);
    }
}
float gridpp::BarnesStructure::localization_distance(float h) const {
    return sqrt(-2*log(m_min_rho)) * h;
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
    float hdist = gridpp::KDTree::calc_straight_distance(p1, p2);
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
    float hdist = gridpp::KDTree::calc_straight_distance(p1, p2);
    if(gridpp::is_valid(m_dist)) {
        if(hdist <= m_dist)
            return 0;
    }
    return m_structure->corr_background(p1, p2);
}
vec gridpp::CrossValidation::corr_background(const Point& p1, const std::vector<Point>& p2) const {
    vec curr_corr = m_structure->corr_background(p1, p2);
    for(int i = 0; i < curr_corr.size(); i++) {
        float hdist = gridpp::KDTree::calc_straight_distance(p1, p2[i]);
        if(gridpp::is_valid(m_dist)) {
            if(hdist <= m_dist)
                curr_corr[i] = 0;
        }
    }
    return curr_corr;
}
gridpp::StructureFunction* gridpp::CrossValidation::clone() const {
    gridpp::StructureFunction* val = new gridpp::CrossValidation(*m_structure);
    return val;
}
float gridpp::CrossValidation::localization_distance(const Point& p) const {
    return m_structure->localization_distance(p);
}
