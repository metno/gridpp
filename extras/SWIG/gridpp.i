%module gridpp
%include "typemaps.i"
%{
#define SWIG_FILE_WITH_INIT
#include <iostream>
void PRINT_DEBUG(std::string message) {
#if 0
    std::cout << message << std::endl;
#endif
}
%}
%include "numpy.i"
%init %{
    import_array();
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
    catch (std::exception &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    }
    catch (...) {
         SWIG_exception(SWIG_RuntimeError, "Unknown exception");
    }
}
%include "std_vector.i"
namespace std {
%template(IntVector) vector<int>;
%template(FloatVector) vector<float>;
%template(FloatVector3) vector<vector<vector<float> > >;
%template(IntVector2) vector<vector<int> >;
%template(DoubleVector) vector<double>;
%template(FloatVector2) vector<vector<float> >;
%template(DoubleVector2) vector<vector<double> >;
}

%apply (float* IN_ARRAY1, int DIM1) {(float* v, int n)}
%apply (double* IN_ARRAY1, int DIM1) {(double* v, int n)}

%define %np_vector_typemaps(DTYPE, NPY_DTYPE)
/*
 * 1D vectors
 */
%typemap(in) std::vector<DTYPE> (std::vector<DTYPE>*ptr, PyArrayObject* py_array){
    PRINT_DEBUG("Typemap(in) std::vector<DTYPE>");
    ptr = NULL;
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        PyObject* py_obj;
        if(array_type($input) == NPY_FLOAT) {
            py_obj = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
        }
        else {
            PyObject* py_obj0 = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_FLOAT), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int size = array_size($input, 0);
        $1 = std::vector<DTYPE>(arg, arg+size);
    }
    else {
        ptr = new std::vector<DTYPE>();
        swig::asptr($input, &ptr);
        if(ptr == NULL) {
            std::cout << "Could not convert type" << std::endl;
            throw std::invalid_argument("Could not convert type");
        }
        $1 = *ptr;
    }
}

%typemap(in) const std::vector<DTYPE> & (std::vector<DTYPE>*ptr, std::vector<DTYPE> temp, PyArrayObject* py_array){
    PRINT_DEBUG("Typemap(in) const std::vector<DTYPE> &");
    ptr = NULL;
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        PyObject* py_obj;
        if(array_type($input) == NPY_FLOAT) {
            py_obj = PyArray_FROMANY($input, NPY_DTYPE, 1, 1, NPY_ARRAY_DEFAULT);
        }
        else {
            PyObject* py_obj0 = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_FLOAT), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int size = array_size($input, 0);
        temp = std::vector<DTYPE>(arg, arg+size);
        $1 = &temp;
    }
    else {
        ptr = new std::vector<DTYPE>();
        swig::asptr($input, &ptr);
        if(ptr == NULL) {
            std::cout << "Could not convert type" << std::endl;
            throw std::invalid_argument("Could not convert type");
        }
        $1 = ptr;
    }
}

%typemap(freearg) std::vector<DTYPE>, const std::vector<DTYPE>&, std::vector<std::vector<DTYPE> >, const std::vector<std::vector<DTYPE> >&, std::vector<std::vector<std::vector<DTYPE> > >, const std::vector<std::vector<std::vector<DTYPE> > >& {
    if(ptr$argnum != NULL) {
        delete ptr$argnum;
    }
}

%typemap(out) std::vector<DTYPE> {
    PRINT_DEBUG("Typemap(out) std::vector<DTYPE>");
    npy_intp dims[1] = {$1.size()};
    $result = PyArray_ZEROS(1, dims, NPY_FLOAT, 0);
    for(long i = 0; i < $1.size(); i++) {
        float* ref = (float*) PyArray_GETPTR1((PyArrayObject*) $result, i);
        ref[0] = $1[i];
    }
}

%typecheck(SWIG_TYPECHECK_INTEGER) std::vector<DTYPE>, const std::vector<DTYPE> & {
    PRINT_DEBUG("typecheck std::vector<DTYPE>");
    if(is_array($input))
        $1 = array_numdims($input) == 1 ? 1 : 0;
    else {
        $1 = 1;
    }
}

/*
 * 2D vectors
 */
%typemap(in) std::vector<std::vector<DTYPE> > (std::vector<std::vector<DTYPE> >*ptr, PyArrayObject* py_array){
    PRINT_DEBUG("Typemap(in) std::vector<std::vector<DTYPE> >");
    ptr = NULL;
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        PyObject* py_obj;
        if(array_type($input) == NPY_FLOAT) {
            py_obj = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
        }
        else {
            PyObject* py_obj0 = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_FLOAT), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        std::vector<std::vector<DTYPE> > temp = std::vector<std::vector<DTYPE> >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i] = std::vector<DTYPE>(arg + i*s1, arg + i*s1 + s1);
        }
        $1 = temp;
    }
    else {
        ptr = new std::vector<std::vector<DTYPE> >();
        swig::asptr($input, &ptr);
        if(ptr == NULL) {
            std::cout << "Could not convert type" << std::endl;
            throw std::invalid_argument("Could not convert type");
        }
        $1 = *ptr;
    }
}

%typemap(in) const std::vector<std::vector<DTYPE> > & (std::vector<std::vector<DTYPE> >*ptr, std::vector<std::vector<DTYPE> > temp, PyArrayObject* py_array){
    PRINT_DEBUG("Typemap(in) const std::vector<std::vector<DTYPE> > &");
    ptr = NULL;
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        PyObject* py_obj;
        if(array_type($input) == NPY_FLOAT) {
            py_obj = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
        }
        else {
            PyObject* py_obj0 = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_FLOAT), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        assert(py_array != NULL);
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        temp = std::vector<std::vector<DTYPE> >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i] = std::vector<DTYPE>(arg + i*s1, arg + i*s1 + s1);
        }
        assert(s1 > 0);
        $1 = &temp;
    }
    else {
        ptr = new std::vector<std::vector<DTYPE> >();
        swig::asptr($input, &ptr);
        if(ptr == NULL) {
            std::cout << "Could not convert type" << std::endl;
            throw std::invalid_argument("Could not convert type");
        }
        $1 = ptr;
    }
}

%typemap(out) std::vector<std::vector<DTYPE> > {
    PRINT_DEBUG("Typemap(out) std::vector<std::vector<DTYPE> >");
    const std::vector<std::vector<DTYPE> >& temp = $1;
    int s0 = temp.size();
    int s1 = 0;
    if(s0 != 0)
        s1 = temp[0].size();
    npy_intp dims[2] = {s0, s1};
    $result = PyArray_ZEROS(2, dims, NPY_FLOAT, 0);
    for(long i = 0; i < s0; i++) {
        for(long j = 0; j < s1; j++) {
            float* ref = (float*) PyArray_GETPTR2((PyArrayObject*) $result, i, j);
            ref[0] = temp[i][j];
        }
    }
}

%typecheck(SWIG_TYPECHECK_INTEGER) std::vector<std::vector<DTYPE> >, const std::vector<std::vector<DTYPE> > & {
    PRINT_DEBUG("typecheck std::vector<std::vector<DTYPE> >");
    if(is_array($input))
        $1 = array_numdims($input) == 2 ? 1 : 0;
    else {
        $1 = 1;
    }
}

/*
 * 3D vectors
 */
%typemap(in) std::vector<std::vector<std::vector<DTYPE> > > (std::vector<std::vector<std::vector<DTYPE> > >*ptr, PyArrayObject* py_array){
    PRINT_DEBUG("Typemap(in) std::vector<std::vector<std::vector<DTYPE> > >");
    ptr = NULL;
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        PyObject* py_obj;
        if(array_type($input) == NPY_FLOAT) {
            py_obj = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
        }
        else {
            PyObject* py_obj0 = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_FLOAT), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        int s2 = array_size($input, 2);
        std::vector<std::vector<std::vector<DTYPE> > > temp = std::vector<std::vector<std::vector<DTYPE> > >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i].resize(s1);
            for(int j = 0; j < s1; j++) {
                temp[i][j] = std::vector<DTYPE>(arg + i * s1 * s2 + j * s2, arg + i * s1 * s2 + (j + 1) * s2);
            }
        }
        $1 = temp;
    }
    else {
        ptr = new std::vector<std::vector<std::vector<DTYPE> > >();
        swig::asptr($input, &ptr);
        if(ptr == NULL) {
            std::cout << "Could not convert type" << std::endl;
            throw std::invalid_argument("Could not convert type");
        }
        $1 = *ptr;
    }
}

%typemap(in) const std::vector<std::vector<std::vector<DTYPE> > > & (std::vector<std::vector<std::vector<DTYPE> > >*ptr, std::vector<std::vector<std::vector<DTYPE> > > temp, PyArrayObject* py_array){
    PRINT_DEBUG("Typemap(in) const std::vector<std::vector<std::vector<DTYPE> > > &");
    ptr = NULL;
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        PyObject* py_obj;
        if(array_type($input) == NPY_FLOAT) {
            py_obj = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
        }
        else {
            PyObject* py_obj0 = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_FLOAT), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        assert(py_array != NULL);
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        int s2 = array_size($input, 2);
        temp = std::vector<std::vector<std::vector<DTYPE> > >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i].resize(s1);
            for(int j = 0; j < s1; j++) {
                temp[i][j] = std::vector<DTYPE>(arg + i * s1 * s2 + j * s2, arg + i * s1 * s2 + (j + 1) * s2);
            }
        }
        assert(s1 > 0);
        $1 = &temp;
    }
    else {
        ptr = new std::vector<std::vector<std::vector<DTYPE> > >();
        swig::asptr($input, &ptr);
        if(ptr == NULL) {
            std::cout << "Could not convert type" << std::endl;
            throw std::invalid_argument("Could not convert type");
        }
        $1 = ptr;
    }
}

%typemap(out) std::vector<std::vector<std::vector<DTYPE> > > {
    PRINT_DEBUG("Typemap(out) std::vector<std::vector<std::vector<DTYPE> > >");
    const std::vector<std::vector<std::vector<DTYPE> > >& temp = $1;
    int s0 = temp.size();
    int s1 = 0;
    if(s0 != 0)
        s1 = temp[0].size();
    int s2 = 0;
    if(s0 != 0 && s1 != 0)
        s2 = temp[0][0].size();
    npy_intp dims[3] = {s0, s1, s2};
    $result = PyArray_ZEROS(3, dims, NPY_FLOAT, 0);
    for(long i = 0; i < s0; i++) {
        for(long j = 0; j < s1; j++) {
            for(long k = 0; k < s2; k++) {
                float* ref = (float*) PyArray_GETPTR3((PyArrayObject*) $result, i, j, k);
                ref[0] = temp[i][j][k];
            }
        }
    }
}

%typecheck(SWIG_TYPECHECK_INTEGER) std::vector<std::vector<std::vector<DTYPE> > >, const std::vector<std::vector<std::vector<DTYPE> > > & {
    PRINT_DEBUG( "typecheck std::vector<std::vector<std::vector<DTYPE> > >");
    if(is_array($input))
        $1 = array_numdims($input) == 3 ? 1 : 0;
    else {
        $1 = 1;
    }
}


//%typecheck(SWIG_TYPECHECK_INTEGER) std::vector<std::vector<DTYPE> > {
//    std::cout << "HERE" << "TYPECHECK 2D" << "DTYPE" << " " << "NPY_DTYPE" << std::endl;
//    $1 = array_numdims($input) == 2 ? 1 : 0;
//}


%enddef
%np_vector_typemaps(int, NPY_INT)
%np_vector_typemaps(long, NPY_LONG)
%np_vector_typemaps(float, NPY_FLOAT)
%np_vector_typemaps(double, NPY_DOUBLE)

%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& output };
%{
#include "gridpp.h"
%}


%include "gridpp.h"
