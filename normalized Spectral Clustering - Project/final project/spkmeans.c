#include <stdlib.h>
#include <stdio.h>
#include "spkmeans.h"

/*
-----------START----------- Assert functions -----------START-----------
*/
void assertNotNullMatrix(double** objTocheck) {
    if(objTocheck == NULL) {
        printf("An Error Has Occured");
        exit(0);
    }
}
void assertNotNullArray(double* objTocheck) {
    if(objTocheck == NULL) {
        printf("An Error Has Occured");
        exit(0);
    }
}

void assertNotNullIntArray(int* objTocheck) {
    if(objTocheck == NULL) {
        printf("An Error Has Occured");
        exit(0);
    }
}
/*
-----------END----------- Assert functions -----------END-----------
*/


/*
-----------START----------- General functions -----------START-----------
*/

double minusZero(double num){
    if (-0.00005 < num && num < 0) {
        num = -0.000001;
    }
    return num;
}
void printMatrix(double **mat, int rowsCnt, int colsCnt) {
    int i, j;
    for (i = 0; i < rowsCnt; i++) {
        for (j = 0; j < colsCnt; j++){
            printf("%.4f", minusZero(mat[i][j]) );
            if (j < colsCnt-1){
                printf(",");
            }
        }
        printf("\n");
    }
}

/*
-----------END----------- General functions -----------END-----------
*/


/*
-----------START----------- allocate and free memory -----------START-----------
*/
/* this function allocate double matrix with the given dimension, return NULL on failure*/
double** allocateMatrix(int rowsCnt, int columnsCnt, int useCalloc) {
    int i;
    double **resultMatrix;
    resultMatrix = useCalloc == 1 ? (double**)calloc(rowsCnt, sizeof(double*)) : (double**)malloc(rowsCnt * sizeof(double*));
    if (resultMatrix == NULL) {
        return NULL;
    }
    for(i = 0 ; i < rowsCnt ; i++) {
        resultMatrix[i] = useCalloc == 1 ? (double*)calloc(columnsCnt, sizeof(double)) : (double*)malloc(columnsCnt * sizeof(double));
        if (resultMatrix[i] == NULL) {
            return NULL;
        }
    }
    return resultMatrix;
}
/*free memory of matrix*/
void freeMatrix(double **matrix, int matDimension) {
    int i;
    if(matrix != NULL) {
        for(i = 0; i < matDimension; i++) {
            free(matrix[i]);
        }
    }
    free(matrix);
}

/*
-----------END----------- allocate and free memory -----------END-----------
*/


/*
-----------START----------- The Weighted Adjacency Matrix -----------START-----------
*/
/* this function get 2 vectors and calculate the euclidean distance of the 2 vectors */
double SQeuclideanDistance(double* v1, double* v2, int vectorDimension){
    int i;
    double norm, tmp;
    double* tempArr;
    norm = 0.0;
    tempArr =  (double *) calloc(vectorDimension,sizeof(double));
    for (i = 0; i < vectorDimension; i++)
        tempArr[i]= v1[i]-v2[i];
    for (i = 0; i < vectorDimension; i++)
    {
        tmp = pow(tempArr[i], 2);
        norm += tmp;
    }
    free(tempArr);
    return norm;
}
/* The Weighted Adjacency Matrix */
double** createWAM(double** vectorsArr, int vectorsCnt, int vectorDimension)
{
    int i, j;
    double** resultMatrix;
    double eucNorm;
    resultMatrix = allocateMatrix(vectorsCnt, vectorsCnt, 0);
    if(resultMatrix == NULL)
    {
        return NULL;
    }

    for(i=0 ; i<vectorsCnt ; i++)
    {
        for(j=i ; j < vectorsCnt ; j++)
        {
            if (i == j) {
                resultMatrix[i][j] = 0;
            }
            else {
                eucNorm = SQeuclideanDistance(vectorsArr[i], vectorsArr[j], vectorDimension);
                resultMatrix[i][j] = exp(-eucNorm / 2.0);
                resultMatrix[j][i] = resultMatrix[i][j];
            }
        }
    }
    return resultMatrix;
}

/*
-----------END----------- The Weighted Adjacency Matrix -----------END-----------
*/

/*
-----------START----------- The Diagonal Degree Matrix -----------START-----------
*/
double** getDDM(double** wam, int matDimension)
{
    int i, j;
    double** resultMatrix;
    double sum;
    resultMatrix = allocateMatrix(matDimension, matDimension, 1);
    if(resultMatrix == NULL){
        return NULL;
    }

    for (i=0 ; i<matDimension ; i++)
    {
        sum = 0.0;
        for (j=0 ; j<matDimension ; j++)
        {
            sum += wam[i][j];
        }
        resultMatrix[i][i] = sum;
    }
    return resultMatrix;
}

/*
-----------END----------- The Diagonal Degree Matrix -----------END-----------
*/


/*
-----------START----------- The Graph Laplacian -----------START-----------
*/

/* The  Graph Laplacian */
double** createGL(double** wam, double** ddm, int matDimension)
{
    int i, j;
    double** resultMatrix;
    resultMatrix = allocateMatrix(matDimension, matDimension, 0);
    if(resultMatrix == NULL) return NULL;
    for (i = 0; i < matDimension; i++) {
        for (j = 0; j < matDimension; j++) {
            resultMatrix[i][j] = ddm[i][j] - wam[i][j];
        }
    }
    return resultMatrix;
}

/*
-----------END----------- The Graph Laplacian -----------END-----------
*/



/*
-----------START----------- The Eigengap Heuristic to determine the number of clusters -----------START-----------
*/

int calcEigengapHeuristic(eigen_values_vectors_obj eigenObj) {
    int i, k;
    double delta, maxDelta;
    maxDelta = -1.0;
    k = 0;
    for(i = 0; i < eigenObj.matDimension/2; i++) {
        delta = fabs(eigenObj.eigenValues[i] - eigenObj.eigenValues[i+1]);
        if(delta > maxDelta) {
            maxDelta = delta;
            k = i + 1;
        }
    }

    return k;
}

/*
-----------END----------- The Eigenap Heuristic to determine the number of clusters -----------END-----------
*/

/*
-----------START----------- Jacobi algorithm -----------START-----------
*/
double** creatACopyOfMatrix(double** mat, int rowsCnt, int colsCnt) {
    int i,j;
    double **result;
    result = allocateMatrix(rowsCnt, colsCnt, 0);
    if(result == NULL) {
        return NULL;
    }
    for(i = 0; i < rowsCnt; i++) {
        for(j = 0; j < colsCnt; j++) {
            result[i][j] = mat[i][j];
        }
    }

    return result;
}
void getPivotIndexes(double** mat, int matDimension, int* Pivot_i, int* Pivot_j) {
    int i,j;
    double pivotValue = -1, tmp;
    for(i = 0; i < matDimension; i++) {
        for(j = matDimension - 1; j >= 0; j--) {
            if (i != j){
                tmp = fabs(mat[i][j]);
                if(tmp > pivotValue) {
                    pivotValue = tmp;
                    *Pivot_i = i;
                    *Pivot_j = j;
                }
            }
        }
    }
}
void findSandC(double** matrix, int i, int j, double *c_val, double *s_val) {
    double s, t, c, theta;
    if(matrix[i][j] == 0){
        theta = 0.0;
        *c_val = 1;
        *s_val = 0;   
        return;
    }
    theta = (matrix[j][j] - matrix[i][i]) / (2 * matrix[i][j]);
    if(theta < 0)
    {
        t = (-1) / (fabs(theta) + pow(pow(theta, 2) + 1.0, 0.5));
    }
    else {
        t = 1 / (fabs(theta) + pow(pow(theta, 2) + 1.0, 0.5));
    }
    c = 1 / pow(pow(t, 2) + 1, 0.5);
    s = t*c;
    *c_val = c;
    *s_val = s;
}
void updateEigenVectors (double** eigenVectors, int matDimension, int i, int j, double c, double s)
{
    int k;
    double ki_val, kj_val;
    for(k = 0; k < matDimension; k++) {
        ki_val = eigenVectors[k][i];
        kj_val = eigenVectors[k][j];
        eigenVectors[k][i] = c*ki_val - s*kj_val;
        eigenVectors[k][j] = s*ki_val + c*kj_val;
    }
}
void updateMatTag(double** mat, double**matTag, int matDimension, int i, int j, double c, double s) {
    int r;
    double newMat_ri, newMat_rj;
    for(r = 0; r < matDimension; r++) {
        if(r != i && r != j) {
            newMat_ri = c*mat[r][i] - s*mat[r][j];
            newMat_rj = c*mat[r][j] + s*mat[r][i];
            matTag[r][i] = newMat_ri;
            matTag[i][r] = newMat_ri;
            matTag[r][j] = newMat_rj;
            matTag[j][r] = newMat_rj;
        }
    }
    matTag[i][j] = 0.0;
    matTag[j][i] = 0.0;
    matTag[i][i] = pow(c,2) * mat[i][i] + pow(s,2) * mat[j][j] - 2.0 * s * c * mat[i][j];
    matTag[j][j] = pow(s,2) * mat[i][i] + pow(c,2) * mat[j][j] + 2.0 * s * c * mat[i][j];
}
double creatOffSquare(double** matrix, int matDimension) {
    int i, j;
    double result = 0.0;
    for (i = 0; i < matDimension; i++) {
        for (j = 0; j < matDimension; j++) {
            if(i != j)
            {
                result += pow(matrix[i][j], 2);
            }
        }
    }
    return result;
}
int checkIfConvergence(double** mat, double** matTag, int matDimension) {
    double eps, offMat, offMatTag;
    eps = 1.0 * pow(10,-5);
    offMat = creatOffSquare(mat, matDimension);
    offMatTag = creatOffSquare(matTag, matDimension);

    if(offMatTag == 0 || (offMat - offMatTag) <= eps ) {
        return 1;
    }
    return 0;
}
double* getEigenValues(double** mat, int matDimension) {
    int i;
    double* eigenValues;

    eigenValues = (double*)calloc(matDimension, sizeof(double));
    if(eigenValues == NULL) {
        return NULL;
    }

    for(i = 0; i < matDimension; i++) {
        eigenValues[i] = mat[i][i];
    }
    return eigenValues;
}
/* return identity matrix of matDimension X matDimension */
double** createIdentityMatrix(int matDimension) {
    double** result;
    int i;
    result = allocateMatrix(matDimension, matDimension, 1);
    if (result == NULL) {
        return NULL;
    }

    for (i = 0; i < matDimension; i++) {
        result[i][i] = 1;
    }
    return result;
}
double** transposeMatrix(double** mat, int rowsCnt, int colsCnt) {
    double** transportedMat;
    int i,j;

    transportedMat = allocateMatrix(rowsCnt, colsCnt, 0);
    if(transportedMat == NULL) {
            return NULL;
    }

    for(i = 0; i < rowsCnt; i++) {
        for(j = 0; j < colsCnt; j++) {
            transportedMat[j][i] = mat[i][j];
        }
    }
    return transportedMat;
}
void inPlaceTranspose(double** mat, int rowsCnt, int colsCnt) {
    double** transportedMat;
    int i, j;
    transportedMat = transposeMatrix(mat, rowsCnt, colsCnt);
    assertNotNullMatrix(transportedMat);
    for(i = 0; i < rowsCnt; i++) {
        for(j = 0; j < colsCnt; j++) {
            mat[i][j] = transportedMat[i][j];
        }
    }
    freeMatrix(transportedMat, rowsCnt);
}
/*get eigen_values_vectors_obj object and create matrix at the 1st row the eigenValues and every other row is an eigenVector (printed as columns) - we do this in order to print it as expected*/
void eigenObjToMatrix(eigen_values_vectors_obj eigenObj, double** eigenMat, int vectorsCnt) {
    int i, j;
    for(i = 0; i < vectorsCnt; i++) {
        eigenMat[0][i] = eigenObj.eigenValues[i];
    }

    for(i = 1; i < vectorsCnt + 1; i++) {
        for(j = 0; j < vectorsCnt; j++) {
            eigenMat[i][j] = eigenObj.eigenVectors[i-1][j];
        }
    }
}

void mergeSortedSection(double* eigenValues, int* sortArr, int start_ind, int med_ind, int end_ind){
    int i, j, k;
    int leftArrLen, rightArrLen;
    double *leftEigenValues, *rightEigenValues;
    int *leftIndexes, *rightIndexes;

    leftArrLen = med_ind - start_ind +1;
    rightArrLen = end_ind - med_ind;

    leftIndexes = (int*)malloc(leftArrLen * sizeof(int));
    rightIndexes = (int*)malloc(rightArrLen * sizeof(int));
    leftEigenValues = (double*)malloc(leftArrLen * sizeof(double));
    rightEigenValues = (double*)malloc(rightArrLen * sizeof(double));

    assertNotNullArray(leftEigenValues);
    assertNotNullArray(rightEigenValues);
    assertNotNullIntArray(leftIndexes);
    assertNotNullIntArray(rightIndexes);

    for(i = 0; i < leftArrLen; i++) {
        leftEigenValues[i] = eigenValues[start_ind + i];
        leftIndexes[i] = sortArr[start_ind + i];
    }

    for(i = 0; i < rightArrLen; i++) {
        rightEigenValues[i] = eigenValues[med_ind + i + 1];
        rightIndexes[i] = sortArr[med_ind + i + 1];
    }

    /*update the original arrays to be sorted between start and end indexes*/
    i = 0, j = 0, k = start_ind;
    while(i < leftArrLen && j < rightArrLen) {
        if(leftEigenValues[i] <= rightEigenValues[j]) {
            eigenValues[k] = leftEigenValues[i];
            sortArr[k] = leftIndexes[i];
            i++;
        }
        else {
            eigenValues[k] = rightEigenValues[j];
            sortArr[k] = rightIndexes[j];
            j++;
        }
        k++;
    }
    while (i < leftArrLen) {
        eigenValues[k] = leftEigenValues[i];
        sortArr[k] = leftIndexes[i];
        i++;
        k++;
    }
    while (j < rightArrLen) {
        eigenValues[k] = rightEigenValues[j];
        sortArr[k] = rightIndexes[j];
        j++;
        k++;
    }

    free(leftIndexes);
    free(leftEigenValues);
    free(rightIndexes);
    free(rightEigenValues);

}

/*recursive function to sort the eigenValues and sortArr - will use the merge sort algorithm*/
void sortEigenValues(double* eigenValues, int* sortArr, int start_ind, int end_ind) {
    int med_ind;
    if(start_ind < end_ind) {
        med_ind = start_ind + ((end_ind - start_ind) / 2);

        sortEigenValues(eigenValues, sortArr, start_ind, med_ind);
        sortEigenValues(eigenValues, sortArr, med_ind + 1, end_ind);

        mergeSortedSection(eigenValues, sortArr, start_ind, med_ind, end_ind);
    }
}

void orderEigenVectorsAccordingToSortedEigenValues(double** eigenVectors, int* sortArr, int vectorsCnt) {
    int i, j;
    double** tmpEigenVectors;
    tmpEigenVectors = allocateMatrix(vectorsCnt, vectorsCnt, 0);
    assertNotNullMatrix(tmpEigenVectors);

    for(i = 0; i < vectorsCnt; i++) {
        for(j = 0; j < vectorsCnt; j++) {
            tmpEigenVectors[i][j] = eigenVectors[sortArr[i]][j];
        }
    }

    /*update the original eigenVectors to be sorted according to the eigenValues sorted*/
    for(i = 0; i < vectorsCnt; i++) {
        for(j = 0; j < vectorsCnt; j++) {
            eigenVectors[i][j] = tmpEigenVectors[i][j];
        }
    }
    freeMatrix(tmpEigenVectors, vectorsCnt);
}

eigen_values_vectors_obj runJacobiAlgorithmUpdated(double** matrix, int matDimension) {
    int i, j, pivot_i,pivot_j, isConvergence, iter;
    double s, c, **mat, **matTag, **eigenVectorsMat;
    eigen_values_vectors_obj eigenObj;
    eigenVectorsMat = createIdentityMatrix(matDimension);
    mat = creatACopyOfMatrix(matrix, matDimension, matDimension);
    matTag = creatACopyOfMatrix(matrix, matDimension, matDimension);
    assertNotNullMatrix(eigenVectorsMat);
    assertNotNullMatrix(mat);
    assertNotNullMatrix(matTag);

    isConvergence = 0;
    iter = 0;
    while(iter < JACOBI_MAX_ITERATIONS && isConvergence == 0) {

        getPivotIndexes(mat, matDimension, &pivot_i, &pivot_j);
        findSandC(mat, pivot_i, pivot_j, &c, &s);
        updateEigenVectors(eigenVectorsMat, matDimension, pivot_i, pivot_j, c, s);
        updateMatTag(mat, matTag, matDimension, pivot_i, pivot_j, c, s);
        isConvergence = checkIfConvergence(mat, matTag, matDimension);
        /*matTag to mat*/
        for (i = 0; i < matDimension; i++)
        {
            for (j = 0; j < matDimension; j++)
                mat[i][j] =  matTag[i][j];

        }
        iter++;
    }

    eigenObj.eigenVectors = eigenVectorsMat;
    eigenObj.matDimension = matDimension;
    eigenObj.eigenValues = getEigenValues(mat, matDimension);
    assertNotNullArray(eigenObj.eigenValues);
    eigenObj.sortArr =  (int*) calloc(matDimension , sizeof (int));
    assertNotNullIntArray(eigenObj.sortArr);

    for (i = 0; i < matDimension; i++) {
        eigenObj.sortArr[i] = i;
    }

    freeMatrix(mat, matDimension);
    freeMatrix(matTag, matDimension);
    return eigenObj;
}

void freeEigenObj(eigen_values_vectors_obj eigenObj){
    free(eigenObj.eigenValues);
    freeMatrix(eigenObj.eigenVectors, eigenObj.matDimension);
    free(eigenObj.sortArr);
    eigenObj.eigenValues = NULL;
    eigenObj.sortArr = NULL;
}

/*
-----------END----------- Jacobi algorithm -----------END-----------
*/


/*
-----------START----------- Reading the input and initialization of the variables -----------START-----------
*/

/* set the goal according to the string from the input */
void setOurGoal(enum goal* g, char *strGoal) {
	if (strcmp(strGoal, "wam") == 0) {
	    *g = wam;
	}
	else if (strcmp(strGoal, "ddg") == 0) {
	    *g = ddg;
	}
	else if (strcmp(strGoal, "gl") == 0) {
	    *g = gl;
	}
	else if (strcmp(strGoal, "jacobi") == 0) {
	    *g = jacobi;
	}
}

/*read line and return length of the line include commas*/
int readDataPoint(FILE *fd, char **line)
{
    int i, dim = 10;
    char * cord, ch;
    cord = (char*) malloc (dim * sizeof(char));
    if(cord == NULL)
    {
        return -1;
    }
    ch = fgetc(fd);
    for(i = 0; ch != '\n' && ch!= EOF ; i++){
        cord[i] = ch;
        if (i + 1 == dim)
        {
            dim *= 2;
            cord = realloc(cord, dim * sizeof(char));
        }
        ch = fgetc(fd);
    }
    cord[i] = '\0';
    cord = realloc(cord, i); /*free not needed memory*/
    *line = cord;
    return i;

}

int getVecDimension(char* line)
{
    int dim = 0, vectorDimension = 0;
    while(line[dim] != '\0')
    {
        if(line[dim] == ',')
        {
            vectorDimension++;
        }
        dim++;
    }
    return ++vectorDimension;
}

double* stringToFloat(char* vecStr, int vectorDimension) {
    char* ch;
    double* vector;
    int counter = 0;
    vector = (double*) malloc (vectorDimension * sizeof(double));
    if(vector == NULL)
    {
        return NULL;
    }
    while(counter < vectorDimension) {
        ch = vecStr + 1;
        while(*ch != ',' && *ch != '\0')
        {
            ch += 1;
        }
        /*replace ',' with '\0' as the end of str */
        *ch = '\0';
        vector[counter] = atof(vecStr);
        counter+=1;
        vecStr = ch + 1;
    }
    return vector;

}

int readInput(FILE *fd, double ***dataPoints, int *vectorDimension) {
    char* line;
    int size, i;
    double* vec;
    *dataPoints = (double**) malloc (1000 * sizeof(double*));
    size = readDataPoint(fd, &line);
    *vectorDimension = getVecDimension(line);
    for(i = 0;line != NULL && line[0] != '\0';i++){
        vec = stringToFloat(line, *vectorDimension);
        if(vec == NULL){
            return -1;
        }
        (*dataPoints)[i] = vec;
        free(line);
        line = NULL;
        size = readDataPoint(fd, &line);
        if(size == -1){
            return -1;
        }
    }
    free(line);

    /*free the memory that not needed*/
    *dataPoints = realloc(*dataPoints, (i) * sizeof(double*));
    return i;
}

/* parse the input and initialize the relevant variables */
int handleInput(int argc, char* argv[], enum goal* g, double *** dataPoints, int *numOfDataPoints, int *vectorDimension) {
    FILE *fd;
    if(argc != 3){
        return 0;
    }

    /*set g (our goal)*/
    setOurGoal(g, argv[1]);

    /*read the file that contain the dataPoints and initialize relevant variable*/
    fd = fopen(argv[2], "r");
    if (fd == NULL) {
        return 0;
    }
    *numOfDataPoints = readInput(fd, dataPoints, vectorDimension);
    if(*numOfDataPoints == -1) {
        return 0;
    }

    fclose(fd);
    return 1;
}

/*
-----------END----------- Reading the input and initialization of the variables -----------END-----------
*/

double** handleWamRequest(double** dataPoints, int numOfDataPoints, int vectorDimension) {
    double **wam;
    wam = createWAM(dataPoints, numOfDataPoints, vectorDimension);
    assertNotNullMatrix(wam);
    return wam;
}

double** handleDDGRequest(double** dataPoints, int numOfDataPoints, int vectorDimension) {
    double **wam, **ddm;
    wam = createWAM(dataPoints, numOfDataPoints, vectorDimension);
    assertNotNullMatrix(wam);

    ddm = getDDM(wam, numOfDataPoints);
    freeMatrix(wam, numOfDataPoints);
    assertNotNullMatrix(ddm);
    return ddm;
}

double** handleGlRequest(double** dataPoints, int numOfDataPoints, int vectorDimension) {
    double **wam, **ddm, **gl;
    wam = createWAM(dataPoints, numOfDataPoints, vectorDimension);
    assertNotNullMatrix(wam);

    ddm = getDDM(wam, numOfDataPoints);
    assertNotNullMatrix(ddm);

    gl = createGL(wam, ddm, numOfDataPoints);
    freeMatrix(wam, numOfDataPoints);
    freeMatrix(ddm, numOfDataPoints);
    assertNotNullMatrix(gl);

    return gl;
}

double** handleJacobiRequest(double** dataPoints, int numOfDataPoints) {
    eigen_values_vectors_obj eigenObj;
    double **jacobi;
    eigenObj = runJacobiAlgorithmUpdated(dataPoints, numOfDataPoints);
    jacobi = allocateMatrix(numOfDataPoints + 1, numOfDataPoints, 0);
    assertNotNullMatrix(jacobi);
    eigenObjToMatrix(eigenObj, jacobi, numOfDataPoints);
    freeEigenObj(eigenObj);
    return jacobi;
}

UAndKObj createUMatrix(double** dataPoints, int numOfDataPoints, int vectorDimension, int k) {
    double **wam, **ddm, **gl;
    double **U;
    int i,j;
    eigen_values_vectors_obj eigenObj;
    UAndKObj u_and_k;

    /*point 1 in the algorithm 1 in the project file, create the weighted adjacency matrix*/
    wam = createWAM(dataPoints, numOfDataPoints, vectorDimension);
    assertNotNullMatrix(wam);

    /*point 2 in the algorithm 1 in the project file, create the diagonal degree matrix*/
    ddm = getDDM(wam, numOfDataPoints);
    assertNotNullMatrix(ddm);

    gl = createGL(wam, ddm, numOfDataPoints);
    freeMatrix(wam, numOfDataPoints);
    freeMatrix(ddm, numOfDataPoints);
    assertNotNullMatrix(gl);

    /*point 3 in the algorithm 1 in the project file, calculate the eigen Values and Vectors of the DDM*/
    eigenObj = runJacobiAlgorithmUpdated(gl, numOfDataPoints);
    freeMatrix(gl, numOfDataPoints);
    inPlaceTranspose(eigenObj.eigenVectors, numOfDataPoints, numOfDataPoints);


    /*sort eigenValues and then sort the eigenVectors*/
    if(k == 0) {
        sortEigenValues(eigenObj.eigenValues, eigenObj.sortArr, 0, numOfDataPoints - 1);
        orderEigenVectorsAccordingToSortedEigenValues(eigenObj.eigenVectors, eigenObj.sortArr, numOfDataPoints);
        k = calcEigengapHeuristic(eigenObj);
    }

    /*point 4 in the algorithm 1 in the project file, each one of the first k eigenVectors will be a column in U*/
    U = allocateMatrix(numOfDataPoints, k, 0);
    assertNotNullMatrix(U);
    for(i = 0; i < numOfDataPoints; i++) {
        for(j = 0; j < k; j++) {
            U[i][j] = eigenObj.eigenVectors[j][i];
        }
    }
    freeEigenObj(eigenObj);
    u_and_k.U = U;
    u_and_k.k = k;

    return u_and_k;
}



int main(int argc, char *argv[])
{
    double **dataPoints = NULL, **wamMat = NULL, **ddmMat = NULL,
     **glMat = NULL, **jacobiMat = NULL;
    enum goal g;
    int numOfDataPoints, vectorDimension, failure = 0;

    if(handleInput(argc, argv, &g, &dataPoints, &numOfDataPoints, &vectorDimension) != 1) {
        printf("An Error Has Occured");
        return 1;
    }
    switch (g)
    {
        case wam:
            wamMat = handleWamRequest(dataPoints, numOfDataPoints, vectorDimension);
            printMatrix(wamMat, numOfDataPoints, numOfDataPoints);
            freeMatrix(wamMat, numOfDataPoints);
            break;
        case ddg:
            ddmMat = handleDDGRequest(dataPoints, numOfDataPoints, vectorDimension);
            printMatrix(ddmMat, numOfDataPoints, numOfDataPoints);
            freeMatrix(ddmMat, numOfDataPoints);
            break;
        case gl:
            glMat = handleGlRequest(dataPoints, numOfDataPoints, vectorDimension);
            printMatrix(glMat, numOfDataPoints, numOfDataPoints);
            freeMatrix(glMat, numOfDataPoints);
            break;
        case jacobi:
            jacobiMat = handleJacobiRequest(dataPoints, numOfDataPoints);
            printMatrix(jacobiMat, numOfDataPoints + 1, numOfDataPoints);
            freeMatrix(jacobiMat, numOfDataPoints);
            break;
        default:
            failure = 1;
            break;
    }

    if(failure == 1) {
        printf("An Error Has Occured\n");
        exit(0);
    }
    freeMatrix(dataPoints, numOfDataPoints);

    return 0;
}
