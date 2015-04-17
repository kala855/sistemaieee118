#include<stdio.h>

typedef struct{
    double *lineas;
    double *gen;
    double *cargas;
    int fnom;
    double numN;
    int numL;
    int numG;
    int numC;
    int maxIter;
} structData;

double maxLineas(structData *data, int widthLineas,int heightLineas);
int calcularYbus(structData *data, double *ybusReal, double * ybusImag);
int ones(int size,double *Vn);
int zeros(int size,double *An);
int printData(structData *data, int widthLineas, int heightLineas, int widthCargas, int heightCargas, int widthGen, int heightGen);
int loadDataFromFile(char *filenameLineas, char *filenameCargas, char* filenameGen, structData *data);
int genVector(int *NNP, int initNumber, int finalNumber);
int setdiff(int *vector1, double *vector2, int size1, int size2, int *c);
int calcularJacobiano(structData *data, double *ybusReal, double *ybusImag, double *Vn, double *An, \
        double *Jpp, double *Jpq, double *Jqp, double *Jqq, double *Pn, double *Qn);
