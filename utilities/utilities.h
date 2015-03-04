#include<stdio.h>

typedef struct{
    double *lineas;
    double *gen;
    double *cargas;
    int numN;
    int numL;
    int numG;
    int numC;
    int maxIter;
} structData;

int loadDataFromFile(char *filename, structData *data);

