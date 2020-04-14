%module gridpp
%include "typemaps.i"
%include "std_vector.i"
%include "std_string.i"
%include exception.i
/*
      SWIG_MemoryError
      SWIG_IOError
      SWIG_RuntimeError
      SWIG_IndexError
      SWIG_TypeError
      SWIG_DivisionByZero
      SWIG_OverflowError
      SWIG_SyntaxError
      SWIG_ValueError
      SWIG_SystemError
*/
%exception {
    try {
        $action
    }
    catch (std::invalid_argument &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_ValueError, s.c_str());
    }
    catch (std::exception &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    }
    catch (...) {
         SWIG_exception(SWIG_RuntimeError, "Unknown exception");
    }
}

%begin %{
#define SWIG_PYTHON_CAST_MODE
    %}
namespace std {
  %template(IntVector) vector<int>;
  %template(IntVector2) vector<vector<int> >;
  %template(DoubleVector) vector<double>;
  %template(FloatVector) vector<float>;
  %template(FloatVector2) vector<vector<float> >;
  %template(FloatVector3) vector<vector<vector<float> > >;
  %template(DoubleVector2) vector<vector<double> >;
}
%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& output };
%{
/*  Put header files here or function declarations like below */
#include "gridpp.h"
%}

%include "gridpp.h"
