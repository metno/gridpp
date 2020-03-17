%module gridpp
%include "typemaps.i"
%include "std_vector.i"
%include "std_string.i"
namespace std {
  %template(IntVector) vector<int>;
  %template(IntVector2) vector<vector<int> >;
  %template(DoubleVector) vector<double>;
  %template(FloatVector) vector<float>;
  %template(FloatVector2) vector<vector<float> >;
  %template(FloatVector3) vector<vector<vector<float> > >;
}
%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& output };
%{
/*  Put header files here or function declarations like below */
#include "gridpp.h"
%}
%include "gridpp.h"
