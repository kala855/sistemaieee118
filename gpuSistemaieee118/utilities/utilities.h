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
int createJacR(int *NNP, int *NNQ, int NumQ, int NumP, int numN, double *Jpp, \
        double *Jpq, double *Jqp, double *Jqq, double *JacR);
int createdPdQ(double *dp, double *dq, int NumP, int NumQ, double *dPdQ);
int transposeJacR(double *JacR,int NumPQ, double *JacRt);
double maxAbs(int NumPQ, double *dPdQ);
int calcCargLineas(structData *data,double *An, double *Vn,double *Ism);
int printDataToFileVec(char *name, int size,double *data);
int printDataToFileMat(char *name, int size,double *data);
__global__ void d_zeros(int size,double *An);
__global__ void d_ones(int size,double *Vn);
__global__ void d_calcularJacobiano_1(int numN, double *ybusReal, double *ybusImag, double *Vn, \
        double *An, double *Pn, double *Qn);
__global__ void d_calcularJacobiano_2(int numN, double *ybusReal, double *ybusImag, double *Vn,\
        double *An,double *Pn, double *Qn, double *Jpp, double *Jpq, double *Jqp, double *Jqq);
__global__ void dp_compute(int NumP,int *NNP, double *Pref, double *Pn, double *dP);
__global__ void dq_compute(int NumQ,int *NNQ, double *Qref, double *Qn, double *dQ);
__global__ void d_createJacR_1(int *NNP, int NumQ, int NumP,int numN, double *Jpp, double *JacR);
__global__ void d_createJacR_2(int *NNP, int *NNQ, int NumQ, int NumP, int numN, double *Jpq, \
        double *JacR);
__global__ void d_createJacR_3(int *NNP, int *NNQ, int NumQ, int NumP, int numN,double *Jpq, \
        double *JacR);
__global__ void d_createJacR_4(int *NNQ, int NumQ, int NumP, int numN, double *Jqq, double *JacR);
__global__ void d_transposeJacr(double *JacR,int NumPQ, double *JacRt);
__global__ void d_filldPdQ(double *d_dQ, int NumQ, int NumP, double *d_dPdQ);
__global__ void d_filldPdQ1(double *d_dP, int NumP, double *d_dPdQ);
__global__ void d_fill_d_dx(double *d_dPdQ, int NumPQ, double *d_dX);
__global__ void d_calc_An(double *dX, int *NNP, int NumP, double *An);
__global__ void d_calc_Vn(double *dX, int *NNQ, int NumP, int NumQ, double *Vn);
