#ifndef PENTMATRIX_H
#define PENTMATRIX_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"

typedef struct {
    PyObject_HEAD
    uint16_t matrix_width;
    uint16_t matrix_height;
    uint32_t modulo_value;
    uint32_t* matrix_values;
} ModuloMatrix;


#endif // PENTMATRIX_H
