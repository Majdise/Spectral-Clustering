#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

PyMODINIT_FUNC PyInit_spkmeans(void);

/*the following 5 functions will be used to move from c to python code and vice versa*/
static PyObject* convert_c_arr_to_python_list(double* c_vector, Py_ssize_t vector_len) {
    Py_ssize_t i;
    PyObject *result;
    PyObject *value;
    result = PyList_New(vector_len);
    for(i = 0; i < vector_len; i++) {
        value = PyFloat_FromDouble(c_vector[i]);
        PyList_SetItem(result,i,value);
    }

    return result;
}

static PyObject* convert_c_matrix_to_python(double** c_matrix, Py_ssize_t rowsCnt, Py_ssize_t colsCnt) {
    Py_ssize_t i;
    PyObject *result;
    PyObject *rowInMatrix;
    result = PyList_New(rowsCnt);
    for(i = 0; i < rowsCnt; i++) {
        rowInMatrix = convert_c_arr_to_python_list(c_matrix[i], colsCnt);
        PyList_SetItem(result,i,rowInMatrix);
    }
    return result;
}

static double* convert_python_lst_to_c_array(PyObject *listObj, Py_ssize_t list_len) {
    Py_ssize_t i;
    PyObject *value;
    double c_value;
    double* result;

    if(!PyList_Check(listObj)) {
        printf("An Error Has Occurred");
        exit(0);
    }
    result = (double*)calloc(list_len, sizeof(double));

    assertNotNullArray(result);

    for(i = 0; i < list_len; i++) {
        value = PyList_GetItem(listObj, i);
        c_value = PyFloat_AsDouble(value);
        if(c_value == -1 && PyErr_Occurred()){
            free(result);
            return NULL;
        }
        result[i] = c_value;
    }

    return result;
}

static double** convert_python_lst_to_c_matrix(PyObject *listObj, Py_ssize_t rowsCnt, Py_ssize_t colsCnt) {

    Py_ssize_t i;
    PyObject *rowInMatrix;
    double **result;

    if(!PyList_Check(listObj)) {
        printf("An Error Has Occurred");
        exit(0);
    }

    result = (double**)calloc(rowsCnt, sizeof(double*));

    assertNotNullMatrix(result);

    for(i = 0; i < rowsCnt; i++) {
        rowInMatrix = PyList_GetItem(listObj, i);
        if(!PyList_Check(rowInMatrix)) {
            printf("An Error Has Occurred");
            exit(0);
        }
        result[i] = convert_python_lst_to_c_array(rowInMatrix, colsCnt);
    }

    return result;
}

/*the following methods will be used to call the C methods*/
static PyObject* wam_matrix(PyObject *self, PyObject *args){
    int data_len, vector_len;
    PyObject *dataPointsPy, *wamMatrixPy;
    double **dataPoints, **wamMatrix;
    if(!PyArg_ParseTuple(args, "O",&dataPointsPy)) {
        printf("An Error Has Occurred");
        return NULL;
    }
    if(!PyList_Check(dataPointsPy)) {
        printf("An Error Has Occurred");
        exit(0);
    }
    data_len = (int) PyList_Size(dataPointsPy);
    vector_len = (int)PyList_Size(PyList_GetItem(dataPointsPy, 0));

    dataPoints = convert_python_lst_to_c_matrix(dataPointsPy, data_len, vector_len);

    wamMatrix = handleWamRequest(dataPoints, data_len, vector_len);

    wamMatrixPy = convert_c_matrix_to_python(wamMatrix, data_len, data_len);

    freeMatrix(dataPoints, data_len);

    return wamMatrixPy;
}

static PyObject* ddg_matrix(PyObject *self, PyObject *args){
    int  data_len, vector_len;
    PyObject *dataPointsPy, *ddgMatrixPy;
    double **dataPoints, **ddgMatrix;
    if(!PyArg_ParseTuple(args, "O", &dataPointsPy)) {
        printf("An Error Has Occurred");
        return NULL;
    }
    if(!PyList_Check(dataPointsPy)) {
        printf("An Error Has Occurred");
        exit(0);
    }
    data_len = (int) PyList_Size(dataPointsPy);
    vector_len = (int)PyList_Size(PyList_GetItem(dataPointsPy, 0));

    dataPoints = convert_python_lst_to_c_matrix(dataPointsPy, data_len, vector_len);

    ddgMatrix = handleDDGRequest(dataPoints, data_len, vector_len);

    ddgMatrixPy = convert_c_matrix_to_python(ddgMatrix, data_len, data_len);

    freeMatrix(dataPoints, data_len);

    return ddgMatrixPy;
}

static PyObject* gl_matrix(PyObject *self, PyObject *args){
    int data_len, vector_len;
    PyObject *dataPointsPy, *glMatrixPy;
    double **dataPoints, **glMatrix;
    if(!PyArg_ParseTuple(args, "O" ,&dataPointsPy)) {
        printf("An Error Has Occurred");
        return NULL;
    }
    if(!PyList_Check(dataPointsPy)) {
        printf("An Error Has Occurred");
        exit(0);
    }
    data_len = (int) PyList_Size(dataPointsPy);
    vector_len = (int)PyList_Size(PyList_GetItem(dataPointsPy, 0));

    dataPoints = convert_python_lst_to_c_matrix(dataPointsPy, data_len, vector_len);

    glMatrix = handleGlRequest(dataPoints, data_len, vector_len);

    glMatrixPy = convert_c_matrix_to_python(glMatrix, data_len, data_len);

    freeMatrix(dataPoints, data_len);

    return glMatrixPy;
}

static PyObject* jacobi_matrix(PyObject *self, PyObject *args) {
    int data_len;
    PyObject *dataPointsPy, *jacobiMatrixPy;
    double **dataPoints, **jacobiMatrix;
    if(!PyArg_ParseTuple(args, "O" ,&dataPointsPy)) {
        printf("An Error Has Occurred");
        return NULL;
    }
    if(!PyList_Check(dataPointsPy)) {
        printf("An Error Has Occurred");
        exit(0);
    }
    data_len = (int) PyList_Size(dataPointsPy);
    dataPoints = convert_python_lst_to_c_matrix(dataPointsPy, data_len, data_len);

    jacobiMatrix = handleJacobiRequest(dataPoints, data_len);

    jacobiMatrixPy = convert_c_matrix_to_python(jacobiMatrix, data_len + 1, data_len);

    freeMatrix(dataPoints, data_len);

    return jacobiMatrixPy;
}

static PyObject* U_matrix(PyObject *self, PyObject *args){
    int k, data_len, vector_len;
    PyObject *dataPointsPy, *UMatrix;
    double **dataPoints;
    UAndKObj u_and_k_obj;
    PyObject* matAndK;
    matAndK = PyList_New( (Py_ssize_t) 2);
    UMatrix = Py_None;
    if(!PyArg_ParseTuple(args, "Oi" ,&dataPointsPy ,&k)) {
        printf("An Error Has Occurred");
                return NULL;
    }
    if(!PyList_Check(dataPointsPy)) {
        printf("An Error Has Occurred");
        exit(0);
    }
    data_len = (int) PyList_Size(dataPointsPy);
    vector_len = (int)PyList_Size(PyList_GetItem(dataPointsPy, 0));

    dataPoints = convert_python_lst_to_c_matrix(dataPointsPy, data_len, vector_len);
    u_and_k_obj = createUMatrix(dataPoints, data_len, vector_len, k);

    UMatrix = convert_c_matrix_to_python(u_and_k_obj.U, data_len, u_and_k_obj.k);

    PyList_SetItem(matAndK,0, UMatrix);
    PyList_SetItem(matAndK,1,Py_BuildValue("i", u_and_k_obj.k));

    freeMatrix(dataPoints, data_len);

    return matAndK;
}


static PyMethodDef spkmeansMethods[] = {
        {"wam", (PyCFunction)wam_matrix,  METH_VARARGS, "receive dataPoints and calculate the Weighted Adjacency Matrix"},
        {"ddg", (PyCFunction)ddg_matrix, METH_VARARGS, "receive dataPoints and calculate the Diagonal Degree Matrix"},
        {"gl", (PyCFunction)gl_matrix, METH_VARARGS, "receive dataPoints and calculate the Diagonal Degree Matrix"},
        {"jacobi", (PyCFunction)jacobi_matrix, METH_VARARGS, "receive dataPoints and calculate the Diagonal Degree Matrix"},
        {"createU", (PyCFunction)U_matrix, METH_VARARGS, "receive dataPoints and calculate the Diagonal Degree Matrix"},

        {NULL, NULL, 0, NULL}
};


static struct PyModuleDef _moduledef = {
         PyModuleDef_HEAD_INIT,
        "spkmeans",
        NULL,
        -1,
        spkmeansMethods
};

PyMODINIT_FUNC PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if(!m) {
        return NULL;
    }
    return m;
}

