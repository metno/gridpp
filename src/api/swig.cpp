#include "gridpp.h"
#include <iostream>

using namespace gridpp;

float* gridpp::test_array(float* v, int n) {
    int count = 0;
    for(int i = 0; i < n; i++)
        count++;
    return v;
 }
float gridpp::test_vec_input(const vec& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++)
        total += input[i];
    return total;
}
int gridpp::test_ivec_input(const ivec& input) {
    int total = 0;
    for(int i = 0; i < input.size(); i++) {
        total += int(input[i]);
    }
    return total;
}
float gridpp::test_vec2_input(const vec2& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++) {
        for(int j = 0; j < input[i].size(); j++) {
            total += input[i][j];
        }
    }
    return total;
}
float gridpp::test_vec3_input(const vec3& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++) {
        for(int j = 0; j < input[i].size(); j++) {
            for(int k = 0; k < input[i][j].size(); k++) {
                total += input[i][j][k];
            }
        }
    }
    return total;
}
vec gridpp::test_vec_output() {
    vec output(3, swig_default_value);
    return output;
}
vec2 gridpp::test_vec2_output() {
    vec2 output(3);
    for(int i = 0; i < 3; i++)
        output[i].resize(3, swig_default_value);
    return output;
}
vec3 gridpp::test_vec3_output() {
    vec3 output(3);
    for(int i = 0; i < 3; i++) {
        output[i].resize(3);
        for(int j = 0; j < 3; j++)
            output[i][j].resize(3, swig_default_value);
    }
    return output;
}

ivec gridpp::test_ivec_output() {
    ivec output(3, swig_default_value);
    return output;
}
ivec2 gridpp::test_ivec2_output() {
    ivec2 output(3);
    for(int i = 0; i < 3; i++)
        output[i].resize(3, swig_default_value);
    return output;
}

ivec3 gridpp::test_ivec3_output() {
    ivec3 output(3);
    for(int i = 0; i < 3; i++) {
        output[i].resize(3);
        for(int j = 0; j < 3; j++)
            output[i][j].resize(3, swig_default_value);
    }
    return output;
}

float gridpp::test_vec_argout(vec& distances) {
    distances.clear();
    distances.resize(10, swig_default_value);
    return 0;
}
float gridpp::test_vec2_argout(vec2& distances) {
    distances.clear();
    distances.resize(10);
    for(int i = 0; i < 10; i++)
        distances[i].resize(10, swig_default_value);
    return 0;
}
void gridpp::test_not_implemented_exception() {
    throw gridpp::not_implemented_exception();
}
vec2 gridpp::test_args_for_R(const gridpp::Points& bpoints,
/*                             const gridpp::StructureFunction& structure,*/
                             const vec2& background) {
    int nS = bpoints.size();
    int nEns = background[0].size();
    int nY = background.size();
    float missing_value = -99999.999;
/*    gridpp::StructureFunction structure = gridpp::BarnesStructure(10000, 0, 0, 0);*/
    BarnesStructure structure = BarnesStructure(10000,1000000,0.5,100000);
    vec2 output = gridpp::init_vec2(nY, nEns, missing_value);
    vec blats = bpoints.get_lats();
    vec blons = bpoints.get_lons();
    vec belevs = bpoints.get_elevs();
    vec blafs = bpoints.get_lafs();
    Point p2 = bpoints.get_point(0);
    for(int i = 0; i < nS; i++) {
        Point p1 = bpoints.get_point(i);
        float localizationRadius = structure.localization_distance(p1);
        float rhos = structure.corr_background(p1, p2);
        std::cout << "----------------------------" << std::endl;
        std::cout << i << " localizationRadius " << localizationRadius << std::endl;
        std::cout << i << " rhos " << rhos << std::endl;
        std::cout << i << " blats " << blats[i] << std::endl;
        std::cout << i << " blons " << blons[i] << std::endl;
        std::cout << i << " belevs " << belevs[i] << std::endl;
        std::cout << i << " blafs " << blafs[i] << std::endl;
        for(int j = 0; j < background[i].size(); j++) {
            std::cout << j << " bkg " << background[i][j] << std::endl;
            output[i][j] = background[i][j];
        }
    }
    return output;
}
