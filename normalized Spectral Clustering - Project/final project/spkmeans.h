#ifndef SPKMEANS_H
#define SPKMEANS_H

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>


#ifndef SPK_CONSTANTS
#define SPK_CONSTANTS
#define JACOBI_MAX_ITERATIONS 100
#endif

#ifndef U_AND_K
#define U_AND_K
typedef struct UAndKObj{
    double **U;
    int k;
} UAndKObj;
#endif


#ifndef EIGEN_VALUES_VECTORS_OBJ
#define EIGEN_VALUES_VECTORS_OBJ
typedef struct eigen_values_vectors_obj {
    int* sortArr;
    double* eigenValues;
    double** eigenVectors;
    int matDimension;
}eigen_values_vectors_obj;
#endif

#ifndef C_ENUM_GOAL
#define C_ENUM_GOAL
enum goal { wam, ddg, gl, jacobi, spk };
#endif


/*
General methods
*/
void freeMatrix(double **matrix, int matDimension);
void assertNotNullArray(double* objTocheck);
void assertNotNullMatrix(double** objTocheck);

/*
Weighted Adjancy Matrix
*/
double** handleWamRequest(double** data, int data_len, int vector_len);
/*
Diagonal Degree Matrix
*/
double** handleDDGRequest(double** data, int data_len, int vector_len);
/*
Graph Laplacian
*/
double** handleGlRequest(double** data, int data_len, int vector_len);
/*
Jacobi algorithm to get eigenvalues and eigenvectors
*/
double** handleJacobiRequest(double** data, int data_len);

/*
services
*/
UAndKObj createUMatrix(double** data, int data_len, int vector_len, int k);


#endif
