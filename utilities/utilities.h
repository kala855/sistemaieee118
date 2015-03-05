#include<stdio.h>

typedef struct{
    double *lineas;
    double *gen;
    double *cargas;
    int fnom;
    int numN;
    int numL;
    int numG;
    int numC;
    int maxIter;
} structData;
int printData(structData *data, int widthLineas, int heightLineas, int widthCargas, int heightCargas, int widthGen, int heightGen);
int loadDataFromFile(char *filenameLineas, char *filenameCargas, char* filenameGen, structData *data);

