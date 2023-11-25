%module gridpp
%init %{
#if defined(SWIGPYTHON)
    import_array();
#endif
    gridpp::initialize_omp();
%}
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
    catch (gridpp::not_implemented_exception &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    }
    catch (std::exception &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    }
    catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown exception");
    }
}
%include "vector.i"
%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& output };
%apply std::vector<std::vector<int> >& OUTPUT { std::vector<std::vector<int> >& count };
%apply std::vector<std::vector<int> >& OUTPUT { std::vector<std::vector<int> >& y_coord };
%apply std::vector<std::vector<int> >& OUTPUT { std::vector<std::vector<int> >& x_coord };
%apply std::vector<float>& OUTPUT { std::vector<float>& distances };
%apply std::vector<float>& OUTPUT { std::vector<float>& standard_error };
%apply std::vector<float>& OUTPUT { std::vector<float>& analysis_variance };
%apply std::vector<float>& OUTPUT { std::vector<float>& output_fcst };
%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& analysis_variance };
%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& distances };
%apply int& OUTPUT { int& X1_out };
%apply int& OUTPUT { int& Y1_out };
%apply int& OUTPUT { int& X2_out };
%apply int& OUTPUT { int& Y2_out };
%apply std::vector<std::vector<std::vector<float> > >& OUTPUT { std::vector<std::vector<std::vector<float> > >& output };
// This turns on automatic description of function arguments in the python package. The
// problem is that it doesn't fetch any of the comments in the header file describing each
// parameter, so only the type is shown.
// A different approach is to use -doxygen with swig
%feature("autodoc", "2");

%{
#include "gridpp.h"
%}
// This is how you add description in the docstring in python. The problem is
// that the comment is added at the end (after the description of function parameters)
// %feature("docstring") gridpp::nearest "Nearest neighbour interpolation"
%include "gridpp.h"
