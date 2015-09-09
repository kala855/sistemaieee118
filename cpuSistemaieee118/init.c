#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utilities/utilities.h"
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
    double *Vn,*An,t;
    structData *data;
    data = (structData*)malloc(sizeof(structData));
    res = loadDataFromFile(fileNameLineas,fileNameCargas,fileNameGen, data);
    data->maxIter = 100;
    data->numN = maxLineas(data,widthLineas,heightLineas);
    data->numL = heightLineas;
    data->numG = heightGen;
    data->numC = heightCargas;
    Vn = (double*)malloc(data->numN*sizeof(double));
    An = (double*)malloc(data->numN*sizeof(double));
    double *ybusReal = (double*) malloc(data->numN*data->numN*sizeof(double));
    double *ybusImag = (double*) malloc(data->numN*data->numN*sizeof(double));
    ones(data->numN,Vn);
    zeros(data->numN,An);
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
    int k;
    int N1;

    for (k = 0; k < data->numG; k++) {
        N1 = (int) data->gen[k*widthGen+0] - 1;
        Pref[N1] = Pref[N1] + data->gen[k*widthGen+1];
        Vn[N1] = data->gen[k*widthGen+2];
    }

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
    int lda = NumP+NumQ,kk;
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

    zeros(data->numL,Ism);
    while (Error>1e-8){
        calcularJacobiano(data,ybusReal,ybusImag,Vn,An,Jpp,Jpq,Jqp,Jqq,Pn,Qn);

        for (i = 0 ; i < NumP ; i++) {
            N1 = NNP[i] - 1;
            dP[i] = Pref[N1] - Pn[N1];
        }

        for (i = 0; i < NumQ; i++ ) {
            N1 = NNQ[i] - 1;
            dQ[i] = Qref[N1] - Qn[N1];
        }

        createJacR(NNP, NNQ, NumQ, NumP, (int)data->numN, Jpp, Jpq, Jqp, Jqq, JacR);
        transposeJacR(JacR,NumPQ,JacRt);
        createdPdQ(dP,dQ,NumP,NumQ,dPdQ);
        memcpy(dX,dPdQ,sizeof(double)*NumPQ);
        dgesv_(&NumPQ,&nrhs,JacRt,&lda,ipiv,dX,&ldb,&info);

        for (k = 0; k < NumP; k++) {
            N1 = NNP[k] - 1;
            An[N1] = An[N1] + dX[k];
        }

        for (k = 0; k < NumQ; k++) {
            N1 = NNQ[k] - 1;
            kk = k + NumP;
            Vn[N1] = Vn[N1] + dX[kk];
        }

       Error = maxAbs(NumPQ,dPdQ);

        if (iter>data->maxIter) {
            printf("..... No converge despues de %d iteraciones\nError = %lf\n", data->maxIter, Error);
            break;
        }
        iter++;
    }



    calcCargLineas(data,An,Vn,Ism);
    printDataToFileVec("ismData",data->numL,Ism);
    printDataToFileVec("vnData",data->numN,Vn);
    printDataToFileVec("anData",data->numN,An);

    free(data);
    free(Vn);
    free(An);
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
