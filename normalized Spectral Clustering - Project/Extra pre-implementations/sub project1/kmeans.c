#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <assert.h>


typedef struct {
    double *lst;
    int l;
} Point;

typedef struct {
    Point centroid;
    double * sumVector;
    int l;
} Cluster;

Point newPoint(char * c);
int countTheD(char * c);
int centroidsFresh(Cluster *c, int k);
void totalVec(Cluster c, Point centro);
double euclidian(Point p1 , Point p2);
Cluster* kmean(int maxIterations, Point* dataPointArr, Cluster* allClusters, int dataPointSize, int k);
void printCentroids(Cluster* finalClusters, int k, int d);
Cluster newClust(Point * pointArr,int x);
int eucNormNoSqrt(double *oldPoints, double *newPoints,Cluster *c,int cnt);
int isNumber(char number[]);

int main(int argc, char *argv[]) {
    char c[1000];
    int max_iter = 200;
    int num, rowLength,j,d_new, k,h,r,rr;
    Point * inputPoints;
    Cluster *allClusters, *finalClusters;
    rowLength = 20;
    inputPoints = (Point *) calloc(rowLength, sizeof(Point));
    num = 0;
    if ((argc<2)){
        printf("An Error Has Occurred");
        return 1;
    }
    if ((argc>3)){
        printf("An Error Has Occurred");
        return 1;
    }
    r = isNumber(argv[1]);
    if (r>0){
        printf("Invalid number of clusters!");
        return 1;
    }
    k = atoi(argv[1]);
    if (k<1){
        printf("Invalid number of clusters!");
        return 1;
    }


    if (argc == 3){
        rr = isNumber(argv[2]);
        if (rr>0){
            printf("Invalid number of iteration!");
            return 1;
        }
        max_iter = atoi(argv[2]);

        if ((max_iter >999)){
            printf("Invalid maximum iteration!");
            return 1;
        }
        if ((max_iter<2)){
            printf("Invalid maximum iteration!");
            return 1;
        }
    }

    if (inputPoints ==NULL){
        printf("An Error Has Occurred");
        return 1;
    } 


   while (fgets(c,1000, stdin) != NULL){
        if (num >= rowLength){
            inputPoints = (Point*) realloc(inputPoints , rowLength * 2 * sizeof(Point));
                if (inputPoints ==NULL){
                printf("An Error Has Occurred");
                return 1;
            }
            rowLength = rowLength * 2;
        }
        inputPoints[num] = newPoint(c);
        num = num + 1;
    }
    if (num<=k){
        printf("Invalid number of clusters!");
        return 1;  
    }


    allClusters = (Cluster *) malloc(k * sizeof(Cluster));
    if (allClusters ==NULL){
        printf("Invalid number of clusters!");
        return 1;
    } 
    for (j = 0 ; j < k; j++) {
        allClusters[j] = newClust(inputPoints, j);
    }
    d_new = inputPoints->l;
    finalClusters = (Cluster *) kmean(max_iter, inputPoints, allClusters, num, k);
    
    
    printCentroids(finalClusters, k, d_new);
    for ( h= 0; h < k; h++){
        free(allClusters[h].sumVector);
        free(allClusters[h].centroid.lst);
        allClusters[h].sumVector = NULL;
        allClusters[h].centroid.lst = NULL;
    }
    for ( h= 0; h < num; h++){
        free(inputPoints[h].lst);
        inputPoints[h].lst= NULL;
    }
    for ( h= 0; h < k; h++){
        free(finalClusters[h].sumVector);
        free(finalClusters[h].centroid.lst);
        finalClusters[h].sumVector = NULL;
        finalClusters[h].centroid.lst = NULL;
    }
    free(inputPoints);
    free(allClusters);
    allClusters = NULL;
    inputPoints = NULL;
    return 0;
}

int isNumber(char number[]){
    int i = 0;
    if (number[0] == '-'){
        i = 1;
    }
    for (; number[i] != 0; i++){
        if (!isdigit(number[i])){
            return 1;
        }
    }
    return 0;
}

Cluster newClust(Point * pointArr,int x){
    Cluster cluster;
    Point point;
    double *lst, *totalVec;
    int e;
    lst = (double *) malloc(pointArr[x].l * sizeof(double));
    for (e = 0 ; e < pointArr[x].l ; e++){
        lst[e] = pointArr[x].lst[e];
    }
    point.lst = lst;
    point.l = pointArr[x].l;
    totalVec = (double *) calloc(pointArr[x].l , sizeof(double));
    cluster.centroid = point;
    cluster.sumVector = totalVec;
    cluster.l = 0;
    return cluster;
}

Point newPoint(char * c){
    int size = countTheD(c) , j = 0;
    const char comma[2] = ",";
    double *float_lst;
    char* tok;
    Point dataPoint;
    float_lst = (double*) calloc(size, sizeof(double));
    tok = strtok(c , comma);
    while (tok != NULL) {
        float_lst[j] = atof(tok);
        tok = strtok(NULL, comma);
        j = j + 1;
    }
    dataPoint.lst = float_lst;
    dataPoint.l = size;
    return dataPoint;
}

int countTheD(char * c){
    int l = 0;
    char *t;
    t = c;
    while (*t != '\0'){
        if (*t == ','){
            l++;
        }
        t++;
    }
    return l+1;
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


Cluster * kmean(int max_iter, Point* inputArr, Cluster* allClusters, int d_new, int k) {
    int countModify = 3;
    int iterations = 0;
    int e,j,closest_Index;
    double min_dis, distance;
    while (iterations < max_iter && countModify>0) {
        iterations= iterations + 1;
        for ( e = 0; e < d_new; e++) {
            min_dis = INT_MAX; 
            closest_Index = -1;
            for ( j = 0 ; j < k; j++) {
                distance = euclidian(allClusters[j].centroid, inputArr[e]);
                if (distance < min_dis) {
                    min_dis = distance;
                    closest_Index = j;
                }
            }
            allClusters[closest_Index].l++;
            totalVec(allClusters[closest_Index], inputArr[e]);
        }
        countModify = centroidsFresh(allClusters, k);
    }
    return allClusters;
}

double euclidian(Point p1 , Point p2){
    double dis;
    int e;
    dis = 0;
    for ( e = 0; e < p1.l; e++){
        dis = dis + ((p1.lst[e] - p2.lst[e])*(p1.lst[e] - p2.lst[e]));
    }
    return dis;
}

void totalVec(Cluster c, Point centro){
    int e;
    for ( e = 0; e < centro.l; e++){
        c.sumVector[e] = c.sumVector[e] + centro.lst[e];
    }
}

int centroidsFresh(Cluster *c, int k){
    double *old_points;
    int cnt=0,e,j,m;
    for (e = 0; e < k; e++){
        old_points = (double*) calloc(c[e].centroid.l, sizeof(double));
        for ( j = 0; j < c[e].centroid.l; j++){
            old_points[j] = c[e].centroid.lst[j];
        }
        if (c[e].l != 0) {
            for (m = 0; m < c[e].centroid.l; m++) {
                c[e].centroid.lst[m] = (c[e].sumVector[m] / (double)c[e].l);
                c[e].sumVector[m] = 0;
            }
        }
        c[e].l = 0;
        cnt = eucNormNoSqrt(old_points, c[e].centroid.lst,c,cnt );
        free(old_points);
    }
    return cnt;  
}

int eucNormNoSqrt(double *oldPoints, double *newPoints,Cluster *c,int cnt){
    double normPow,epsilonPowerTwo;
    int e;
    normPow = 0;
    epsilonPowerTwo= 0.001*0.001;
    for ( e=0; e< c[e].centroid.l;e++){
        normPow = normPow +((newPoints[e])-(oldPoints[e]))*((newPoints[e])-(oldPoints[e]));
    }
    if (normPow >= epsilonPowerTwo){
        cnt = cnt+1;
    }
    return cnt;
}

