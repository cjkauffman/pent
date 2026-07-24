#include "pentmatrix.h"
#include<inttypes.h>
#include <stdio.h>

static PyObject *ModuloMatrix_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    ModuloMatrix *self;
    self = (ModuloMatrix*) type->tp_alloc(type, 0);
    if(self == NULL) {
    //    PyErr_SetString(PyExc_AttributeError, "Modulo Matrix.");
        return NULL;
    }
    self->matrix_width = 0;
    self->matrix_height = 0;
    self->modulo_value = 0;
    self->matrix_values = NULL;
    return (PyObject *) self;
}

static int ModuloMatrix_init(PyObject *op, PyObject *args, PyObject *kwds) {
    ModuloMatrix *self = (ModuloMatrix *) op;
    PyObject *basearray = NULL;
    if (!PyArg_ParseTuple(args, "OI", &basearray, &self->modulo_value)) {return -1;}
    if (!PyList_Check(basearray)) {return -2;}
    self->matrix_height = PyList_Size(basearray);
    if(self->matrix_height>0) {
        if(!PyList_Check(PyList_GetItem(basearray, 0))) {
            return -3;
        } else {
            self->matrix_width = PyList_Size(PyList_GetItem(basearray, 0));
        }
        for (int i = 1; i < self->matrix_height; i++) {
            PyObject* row = PyList_GetItem(basearray, i);
            if(!PyList_Check(row)) {
                return -3;
            }
            if (PyList_Size(row) != self->matrix_width) {
                return -4;
            }
        }
        self->matrix_values = malloc(sizeof(uint32_t)*self->matrix_width*self->matrix_height);
        if (!self->matrix_values){return -5;}
        for (int i = 0; i < self->matrix_height; i++) {
            PyObject* row = PyList_GetItem(basearray, i);
            for (int j = 0; j < self->matrix_width; j++) {
                PyObject* pyvalue = PyList_GetItem(row, j);
                if (!PyLong_Check(pyvalue)) {
                    free(self->matrix_values);
                    return -5;
                }
                int64_t *signedvalue;
                if(PyLong_AsUInt64(pyvalue, signedvalue)) {
                    free(self->matrix_values);
                    return -6;
                }
                *signedvalue %= self->modulo_value;
                if(*signedvalue < 0) {*signedvalue += self->modulo_value;}
                uint32_t matrix_value = *signedvalue & 0xFFFFFFFF;
                self->matrix_values[i*self->matrix_width + j] = matrix_value;
            }

        }
    }
    return 0;
}

static void ModuloMatrix_dealloc(PyObject *op) {
    ModuloMatrix *self = (ModuloMatrix *)op;
    free(self->matrix_values);
    Py_TYPE(self)->tp_free(self);
}


static PyObject *ModuloMatrix_getarray(PyObject *op, PyObject *Py_UNUSED(dummy)) {
    ModuloMatrix *self = (ModuloMatrix *)op;
    PyObject *rowlist = PyList_New(self->matrix_height);
    for(int i = 0; i < self->matrix_height; i++) {
        PyObject *row = PyList_New(self->matrix_width);
        PyList_SetItem(rowlist, i, row);
        for(int j = 0; j < self->matrix_width; j++) {
            PyObject *entry = PyLong_FromUInt32(self->matrix_values[i*self->matrix_width + j]);
            PyList_SetItem(row, j, entry);
        }
    }
    return rowlist;
}

static PyMethodDef ModuloMatrix_methods[] = {
    {"array", ModuloMatrix_getarray, METH_NOARGS, "Returns the associated array."},
    {NULL}
};

static PyTypeObject PyModMatrix = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pent.ModuloMatrix",
    .tp_doc = PyDoc_STR("Modulo Matrix"),
    .tp_basicsize = sizeof(ModuloMatrix),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = ModuloMatrix_new,
    .tp_init = ModuloMatrix_init,
    .tp_dealloc = ModuloMatrix_dealloc,
    .tp_methods = ModuloMatrix_methods
};

static PyObject* ModuloMatrix_product(PyObject* self, PyObject* args){
    PyObject *A, *B;
    if (!PyArg_ParseTuple(args, "OO", &A, &B)) {return NULL;}
    ModuloMatrix *a = (ModuloMatrix *)A;
    ModuloMatrix *b = (ModuloMatrix *)B;
    if(a->matrix_width != b->matrix_height || a->modulo_value != b->modulo_value) {return NULL;}
    ModuloMatrix* product = (ModuloMatrix *)ModuloMatrix_new(&PyModMatrix, NULL, NULL);
    product->matrix_height = a->matrix_height;
    product->matrix_width = b->matrix_width;
    product->modulo_value = a->modulo_value;
    product->matrix_values = malloc(a->matrix_height * b->matrix_width * sizeof(uint32_t));
    for(uint16_t i = 0; i < product->matrix_height; i++) {
        for (uint16_t j = 0; j < product->matrix_width; j++) {
            uint64_t running_sum = 0;
            for(int k = 0; k < a->matrix_width; k++) {
                running_sum += (uint64_t)(a->matrix_values[i*a->matrix_width + k]) * (uint64_t)(b->matrix_values[k*b->matrix_width + j]);
                running_sum %= product->modulo_value;
            }
            product->matrix_values[i*product->matrix_width + j] = (uint32_t)(running_sum & 0xFFFFFFFF);
        }
    }
    return (PyObject *)product;
}



static PyObject* ModuloMatrix_power(PyObject* self, PyObject* args){
    PyObject *A;
    unsigned long long *power;
    if (!PyArg_ParseTuple(args, "OK", &A, power)) {return NULL;}
    ModuloMatrix *a = (ModuloMatrix *)A;
    if(a->matrix_width != a->matrix_height){return NULL;}
    ModuloMatrix* runningproduct = (ModuloMatrix *)ModuloMatrix_new(&PyModMatrix, NULL, NULL);
    runningproduct->matrix_height = runningproduct->matrix_width = a->matrix_width;
    runningproduct->modulo_value = a->modulo_value;
    runningproduct->matrix_values = malloc(a->matrix_height * a->matrix_width * sizeof(uint32_t));
    for(uint32_t i = 0; i < a->matrix_height*a->matrix_width; i++) {
        if(i % (a->matrix_height + 1)){runningproduct->matrix_values[i] = 0;} else {runningproduct->matrix_values[i] = 1;}
    }
    printf("Check 1, power = %llu\n", *power);
    ModuloMatrix* binarypower = (ModuloMatrix *)ModuloMatrix_new(&PyModMatrix, NULL, NULL);
    binarypower->matrix_height = binarypower->matrix_width = a->matrix_height;
    binarypower->modulo_value = a->modulo_value;
    binarypower->matrix_values = malloc(a->matrix_height * a->matrix_width * sizeof(uint32_t));
    for(uint32_t i = 0; i < a->matrix_height*a->matrix_width; i++) {
        binarypower->matrix_values[i] = a->matrix_values[i];
    }
    printf("Check 2, power = %llu\n", *power);
    for(unsigned long long powermask = 1;powermask;powermask *= 2){
        printf("Check 3\n");
        if (powermask & *power) {
            printf("Check 4\n");
            PyObject* rpargs = PyTuple_New(2);
            printf("Check 5\n");
            PyTuple_SetItem(rpargs, 0, (PyObject * )runningproduct);
            printf("Check 6\n");
            PyTuple_SetItem(rpargs, 1, (PyObject * )binarypower);
            printf("Check 7\n");
            ModuloMatrix* newrunningproduct = (ModuloMatrix *)ModuloMatrix_product(self, rpargs);
            runningproduct = newrunningproduct;
        }
        printf("Check 8\n");
        PyObject* bpargs = PyTuple_New(2);
        PyTuple_SetItem(bpargs, 0, (PyObject * )binarypower);
        PyTuple_SetItem(bpargs, 1, (PyObject * )binarypower);
        ModuloMatrix* newbinarypower = (ModuloMatrix *)ModuloMatrix_product(self, bpargs);
        binarypower = newbinarypower;
    }
    return (PyObject *)runningproduct;
}

static PyMethodDef pent_methods[] = {
    {"product", ModuloMatrix_product, METH_VARARGS, "Returns the product of two matrices."},
    {"power", ModuloMatrix_power, METH_VARARGS, "Returns a matrix to the power of a nonnegative integer."},
    {NULL}
};

static int pent_module_exec(PyObject *m) {
    if(PyType_Ready(&PyModMatrix) < 0) {return -1;}
    if (PyModule_AddObjectRef(m, "ModuloMatrix", (PyObject*) &PyModMatrix) < 0) {return -1;}
    return 0;
}

static PyModuleDef_Slot pent_module_slots[] = {
    {Py_mod_exec, pent_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_MULTIPLE_INTERPRETERS_NOT_SUPPORTED},
    {0, NULL}
};
static PyModuleDef pent_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "pent",
    .m_doc = "",
    .m_size = 0,
    .m_slots = pent_module_slots,
    .m_methods = pent_methods,
};
PyMODINIT_FUNC PyInit_pent(void) {
    return PyModuleDef_Init(&pent_module);
}

void printmatrix(ModuloMatrix* a) {
    for(int i = 0; i < a->matrix_height; i++) {
        for(int j = 0; j < a->matrix_width; j++) {
            printf("%"PRIu32" ", a->matrix_values[i*a->matrix_width + j]);
        }
        printf("\n");
    }
}