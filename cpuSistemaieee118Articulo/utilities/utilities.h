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
    int Sbase;
    int iteraciones;
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
int createJacR(int *NNP, int *NNQ, int NumQ, int NumP, int numN, double *Jpp, \
        double *Jpq, double *Jqp, double *Jqq, double *JacR);
int createdPdQ(double *dp, double *dq, int NumP, int NumQ, double *dPdQ);
int transposeJacR(double *JacR,int NumPQ, double *JacRt);
double maxAbs(int NumPQ, double *dPdQ);
int calcCargLineas(structData *data,double *An, double *Vn,double *Ism);
int printDataToFileVec(char *name, int size,double *data);
int printDataToFileMat(char *name, int size,double *data);
int loadCorrientesMax(char *fileNameIMax, double *Imax);
int calcularMatrizA(structData *data, int widthLineas, double *A);
