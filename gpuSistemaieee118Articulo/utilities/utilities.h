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


typedef struct{
    double *sum;
    double *sumcuad;
    double *lv;
    double *hv;
    double *sob;
} Mont;

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
int printMatrixToFile(double *A, int numFilas, int numColumnas, char *name);
int calcularZp(structData *data, int heightLineas, int widthLineas, double *ZpReal,double *ZpImag);
int loadNW(char *fileNameNW, double *NW);
int calculoIrama(double *Vrama, double *ZpReal, double *ZpImag, int heightLineas, double *Irama);
int newtonRaphson(structData *data, double *Vn, double *An, double *ybusReal, double *ybusImag);
int calculoVn(double *Vn, double *An, int height, double *VnReal, double *VnImag);
int addVectors(double *vec1, double *vec2, int numN);
int subVectors(double *vec1, double *vec2, int numN);
int calculoSobrecarga(double *IlineaReal, double *IlineaImag, double *sobrecarga, double *Imax, int numL);
int calculoMontSobrecarga(int numL, double *sobrecarga, Mont *mont);
int calculoCorrientesRama(int numN,double *Vn, Mont *mont);
int calculosFinales(int numN,int ni, Mont *mont, double *Vmedia, double *Vdesv, double *Probmin, double *ProbMax);
int calculoProbSobrecarga(int numL, int ni, Mont *mont, double *Probsobrecarga);
