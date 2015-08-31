#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utilities/utilities.h"
#include <time.h>
#include <math.h>
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,\
        double* b, int* ldb, int* info );


int main(){
    int res,i,j;
    int widthCargas = 3;
    int widthGen = 3;
    int widthLineas = 6;
    int heightCargas = 83;
    int heightGen = 15;
    int heightLineas = 186, NumP;
    char *fileNameLineas = "../../inputs/lineas";
    char *fileNameCargas = "../../inputs/cargas";
    char *fileNameGen = "../../inputs/gen";
    char *fileNameIMax = "../../inputs/imax";
    char *fileNameNW = "../../inputs/NW";
    int numDataImax = 186, NumW = 3;
    double *Imax, *A, *ZpReal, *ZpImag, *NW;
    structData *data;
    Imax = malloc(numDataImax*sizeof(double));
    data = (structData*)malloc(sizeof(structData));
    ///////// Los datos son cargados desde los archivos de texto //////////////////////////
    //////// Lineas, generadores y cargas se llevan a la estructura data//////////////////
    res = loadDataFromFile(fileNameLineas,fileNameCargas,fileNameGen, data);
    data->Sbase = 100;
    data->iteraciones = 40;
    data->fnom = 60.0;
    data->maxIter = 100;
    data->numN = maxLineas(data,widthLineas,heightLineas);
    data->numL = heightLineas;
    data->numG = heightGen;
    data->numC = heightCargas;

    //////////////////////////////////////////////////////////////////////////////////////
    A = malloc((int)(data->numN)*(int)(data->numL)*sizeof(double));
    NW = malloc(3*4*sizeof(double));
    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////
    res = loadCorrientesMax(fileNameIMax, Imax);
    res = loadNW(fileNameNW, NW);

/*    for (i = 0; i < 3; i++) {
       for (j = 0; j < 4; j++) {
           printf("%lf ", NW[i*4+j]);
       }
       printf("\n");
    }*/

    //////////////////////////////////////////////////////////////////////////////////////
    calcularMatrizA(data,widthLineas,A);
    printMatrixToFile(A,data->numL, data->numN, "A");
    //////////////////////////////////////////////////////////////////////////////////////
    ZpReal = malloc(heightLineas*sizeof(double));
    ZpImag = malloc(heightLineas*sizeof(double));
    calcularZp(data,heightLineas,widthLineas,ZpReal,ZpImag);

    int ni = 1, m;
    int k;
    double r, lambda = 0.0, kk, N, Pmax, Vw, Pw;
    NumW = 3;
    double Vmin = 4, Vnom = 12, Vmax = 25;
    srand((unsigned)time(NULL));
    for (k = 0; k < ni; k++) {
        for (m = 0; m < NumW; m++) {
            r = ((double)rand()/(double)RAND_MAX);
            lambda = NW[m*NumW+2];
            kk = NW[m*NumW+3];
            Pmax = NW[m*NumW+1];
            N = NW[m*NumW+0];
            Vw = lambda*pow((-log(1-r)),(1/kk));
            if ((Vw<Vmin)||(Vw>Vmax))
                Pw = 0;
            if ((Vmin<Vw)&&(Vw<Vnom))
                Pw = Pmax*pow((Vw/Vnom),3);
            if ((Vw>Vnom)&&(Vw<Vmax))
                Pw = Pmax;
            printf("%lf \n", Vw);
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////

    double *ybusReal = (double*) malloc(data->numN*data->numN*sizeof(double));
    double *ybusImag = (double*) malloc(data->numN*data->numN*sizeof(double));
    calcularYbus(data,ybusReal,ybusImag);
    NumP = (int) data->numN - 1;
    int *NNP = (int *)malloc((data->numN-1)*sizeof(int));
    genVector(NNP, 2,data->numN);
    int *vector1 = (int*)malloc(data->numN*sizeof(int));
    genVector(vector1, 1, data->numN);
    int *NNQ = (int *) malloc(data->numN*sizeof(int));
    int NumQ = setdiff(vector1, data->gen, data->numN, data->numG, NNQ);
    double *Pref = (double*)malloc(data->numN*sizeof(double));
    double *Qref = (double*)malloc(data->numN*sizeof(double));
    zeros(data->numN,Pref);
    zeros(data->numN,Qref);
    int N1;

    for (k = 0 ; k < data->numC; k++) {
        N1 = (int)data->cargas[k*widthCargas] - 1;
        Pref[N1] = Pref[N1] - data->cargas[k*widthCargas+1];
        Qref[N1] = Qref[N1] - data->cargas[k*widthCargas+2];
    }


    double *dP = (double*)malloc(NumP*sizeof(double));
    double *dQ = (double*)malloc(NumQ*sizeof(double));

    zeros(NumP,dP);
    zeros(NumQ,dQ);

    double Error = 100.0;
    int iter = 0;
    int lda = NumP+NumQ;
    int NumPQ = NumP+NumQ, nrhs = 1;
    int *ipiv,ldb = NumQ+NumP,info;

    double *Jpp, *Jpq, *Jqp, *Jqq, *Pn, *Qn, *JacR, *dPdQ, *JacRt,*dX, *Ism;

    Jpp = (double*)malloc(data->numN*data->numN*sizeof(double));
    Jpq = (double*)malloc(data->numN*data->numN*sizeof(double));
    Jqp = (double*)malloc(data->numN*data->numN*sizeof(double));
    Jqq = (double*)malloc(data->numN*data->numN*sizeof(double));
    Pn = (double*)malloc(data->numN*sizeof(double));
    Qn = (double*)malloc(data->numN*sizeof(double));
    JacR = (double*)malloc((NumPQ)*(NumPQ)*sizeof(double));
    dPdQ = (double*)malloc((NumPQ)*sizeof(double));
    dX = (double*)malloc((NumPQ)*sizeof(double));
    ipiv = (int*)malloc((NumP+NumQ)*sizeof(int));
    JacRt = (double*)malloc((NumPQ)*(NumPQ)*sizeof(double));
    Ism = (double*)malloc(data->numL*sizeof(double));

    printDataToFileVec("ismData",data->numL,Ism);

    free(data);
    free(ybusReal);
    free(ybusImag);
    free(NNP);
    free(NNQ);
    free(vector1);
    free(dP);
    free(dQ);
    free(Pref);
    free(Qref);
    free(JacR);
    free(ipiv);
    free(JacRt);
    free(dX);
    free(Jpp);
    free(Jpq);
    free(Jqp);
    free(Jqq);
    free(Ism);
    return res;
}
