#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

typedef struct {
    double *lst;
    int d;
} Point;

typedef struct {
    Point centroid;
    double * sumVector;
    int d;
} Cluster;

Point create_Point(double *pointVec, int d);
int updateCent(Cluster *c, int k);
void sumVec(Cluster c, Point p);
double eucDis(Point p1 , Point p2);
Cluster* kmean(int maxIterations, Point* dataPointArr, Cluster* clusters, int dataPointSize, int k);
void printCentroids(Cluster* finalClusters, int k, int d);
void printPoints(Point* inputPoint, int num_of_lines, int d);
Cluster initCluster(double *pointVec, int d);
int centFun(int i, int cnt, double *prev_cent, Cluster *c );

static PyObject *kmeans(int, int, int, int, PyObject *, PyObject *);

static PyObject *fit(PyObject *self, PyObject *args);

#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef kmeansMethods[] = {
   FUNC(METH_VARARGS, fit, "calculate k-means"),
           {NULL, NULL, 0, NULL}   /* sentinel */
};
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "mykmeanssp", NULL, -1, kmeansMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject* pyObj;
    pyObj = PyModule_Create(&moduledef);
    if (!pyObj)  {
        return NULL;
    }
    return pyObj;
};

static PyObject *fit(PyObject *self, PyObject *args) {
    int k;
    int max_iter;
    int num_of_lines;
    int dim;
    PyObject *centroids_pyObj;
    PyObject *points_to_cluster_pyObj;

        if (!PyArg_ParseTuple(args, "llllOO:fit", &k, &max_iter, &num_of_lines, &dim, &centroids_pyObj,
                              &points_to_cluster_pyObj)) {
            return NULL;
        }

    return Py_BuildValue("O", kmeans(k, max_iter, num_of_lines, dim, centroids_pyObj, points_to_cluster_pyObj));
}

static PyObject *kmeans(int k, int max_iter, int num_of_lines, int vector_len, PyObject *centroids_pyObj,
                        PyObject *points_to_cluster_pyObj) {
        int i, j, centroids_length, points_to_cluster_length;
        double *points_to_cluster;
        Point *inputPoints;
        double *centroids;
        Cluster *clusters, *final_clusters;

        centroids_length = k * vector_len;
        points_to_cluster_length = num_of_lines * vector_len;

        points_to_cluster = (double *) calloc(points_to_cluster_length, sizeof(double));
        assert(points_to_cluster != NULL && "Can't allocate memory for points_to_cluster, please try again!");
        inputPoints = (Point *) calloc(num_of_lines, sizeof(Point));
        assert(inputPoints != NULL && "Can't allocate memory for inputPoints, please try again!");

        centroids = (double *) calloc(centroids_length, sizeof(double));
        assert(centroids != NULL && "Can't allocate memory for centroids, please try again!");
        clusters = (Cluster *) calloc(k, sizeof(Cluster));
        assert(clusters != NULL && "Can't allocate memory for clusters, please try again!");



        for (i = 0; i < points_to_cluster_length; i++) {
            PyObject *item;
            item = PyList_GetItem(points_to_cluster_pyObj, i);
            points_to_cluster[i] = PyFloat_AsDouble(item);
        }

        for (i = 0; i < num_of_lines; i++) {
            double *pointVec = (double*) calloc(vector_len, sizeof(double));
            for (j = 0; j < vector_len; j++) {
                pointVec[j] = points_to_cluster[i * vector_len + j];
            }

            inputPoints[i] = create_Point(pointVec, vector_len);
        }

        for (i = 0; i < centroids_length; i++) {
            PyObject *item;
            item = PyList_GetItem(centroids_pyObj, i);
            centroids[i] = PyFloat_AsDouble(item);
        }

        for (i = 0; i < k; i++) {
            double *pointVec = (double*) calloc(vector_len, sizeof(double));
            for(j = 0; j < vector_len ; j++) {
                pointVec[j] = centroids[i * vector_len + j];
            }

            clusters[i] = initCluster(pointVec, vector_len);
        }

        final_clusters = (Cluster *) kmean(max_iter, inputPoints, clusters, num_of_lines, k);

        PyObject *list = PyList_New(k * vector_len);
        for (i = 0 ; i< k ; i++) {
            for(j = 0 ; j < vector_len; j++){
                PyList_SetItem(list, (i * vector_len + j) , PyFloat_FromDouble(final_clusters[i].centroid.lst[j]));
            }
        }

        // free memory
        for (i = 0; i < k; i++){
            free(clusters[i].sumVector);
            free(clusters[i].centroid.lst);
            clusters[i].sumVector = NULL;
            clusters[i].centroid.lst = NULL;
        }
        for (i = 0; i < num_of_lines; i++){
            free(inputPoints[i].lst);
            inputPoints[i].lst= NULL;
        }
        for (i = 0; i < k; i++){
            free(final_clusters[i].sumVector);
            free(final_clusters[i].centroid.lst);
            final_clusters[i].sumVector = NULL;
            final_clusters[i].centroid.lst = NULL;
        }
        free(points_to_cluster);
        free(centroids);
        free(inputPoints);
        free(clusters);
        points_to_cluster = NULL;
        centroids = NULL;
        inputPoints = NULL;
        clusters = NULL;
        return list;
}


Cluster initCluster(double *pointVec, int d) {
    Cluster cluster;
    Point point;
    double *lst, *sumVec;
    int i;
    lst = (double *) malloc(d * sizeof(double));
    for (i = 0 ; i < d ; i++){
        lst[i] = pointVec[i];
    }
    point.lst = lst;
    point.d = d;
    sumVec = (double *) calloc(d , sizeof(double));
    cluster.centroid = point;
    cluster.sumVector = sumVec;
    cluster.d = 0;
    return cluster;
}

Point create_Point(double *pointVec, int d) {
    Point point1;
    point1.lst = pointVec;
    point1.d = d;
    return point1;
}

void printPoints(Point* inputPoint, int num_of_lines, int d) {
    int a,b;
    for ( a = 0 ; a < num_of_lines ; a++) {
        for ( b = 0 ; b < d  ; b++) {
            printf("%.4f", inputPoint[a].lst[b]);
            if (b != d-1){
                printf(",");
            }
        }
        printf("\n");
    }
}


void printCentroids(Cluster* finalClusters, int k, int d) {
    int a,b;
    for ( a = 0 ; a < k ; a++) {
        for ( b = 0 ; b < d  ; b++) {
            printf("%.4f", finalClusters[a].centroid.lst[b]);
            if (b != d-1){
                printf(",");
            }
        }
        printf("\n");
    }
}

Cluster * kmean(int max_iter, Point* inputArr, Cluster* clusters, int d_size, int k) {
    int centroidsChanged = 3, iterations = 0;
    int i,j,closest_Index;
    double min_dis, distance;
    while (iterations < max_iter && centroidsChanged > 0 ) {
        iterations= iterations + 1;
        for ( i = 0; i < d_size; i++) {
            min_dis = INT_MAX;
            closest_Index = -1;
            for ( j = 0 ; j < k; j++) {
                distance = eucDis(clusters[j].centroid, inputArr[i]);
                if (distance < min_dis) {
                    min_dis = distance;
                    closest_Index = j;
                }
            }
            clusters[closest_Index].d++;
            sumVec(clusters[closest_Index], inputArr[i]);
        }
        centroidsChanged = updateCent(clusters, k);
    }
    return clusters;
}

double eucDis(Point p1 , Point p2){
    double dis;
    int i;
    dis = 0;
    for ( i = 0; i < p1.d; i++){
        dis = dis + ((p1.lst[i] - p2.lst[i])*(p1.lst[i] - p2.lst[i]));
    }
    return dis;
}

void sumVec(Cluster c, Point p){
    int i;
    for ( i = 0; i < p.d; i++){
        c.sumVector[i] = c.sumVector[i] + p.lst[i];
    }
}

int updateCent(Cluster *c, int k){
    int cnt=0,i,j,m;
    double *prev_cent;
    for (i = 0; i < k; i++){
        prev_cent = (double*) calloc(c[i].centroid.d, sizeof(double));
        for ( j = 0; j < c[i].centroid.d; j++){
            prev_cent[j] = c[i].centroid.lst[j];
        }
        if (c[i].d != 0) {
            for (m = 0; m < c[i].centroid.d; m++) {
                c[i].centroid.lst[m] = (c[i].sumVector[m] / (double)c[i].d);
                c[i].sumVector[m] = 0;
            }
        }
        c[i].d = 0;
        cnt = centFun(i,cnt,prev_cent,c);
        free(prev_cent);
    }
    return cnt;
}
int centFun(int i, int cnt, double *prev_cent, Cluster *c ){
    int n;
    for (n = 0; n < c[i].centroid.d; n++) {
            if (prev_cent[n] != c[i].centroid.lst[n]) {
                cnt = cnt + 1;
                break;
            }
        }
    return cnt;
}
