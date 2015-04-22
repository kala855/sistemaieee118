#include <stdio.h>
#include <stdlib.h>
#include "utilities/utilities.h"

int main(){
    int res;
    int widthCargas = 3;
    int widthGen = 3;
    int widthLineas = 6;
    int heightCargas = 83;
    int heightGen = 15;
    int heightLineas = 186, NumP;
    char *fileNameLineas = "../inputs/lineas";
    char *fileNameCargas = "../inputs/cargas";
    char *fileNameGen = "../inputs/gen";
    double *Vn,*An;
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

    int Error = 100;
    int iter = 0;

    double *Jpp, *Jpq, *Jqp, *Jqq, *Pn, *Qn, *JacR;

    Jpp = (double*)malloc(data->numN*data->numN*sizeof(double));
    Jpq = (double*)malloc(data->numN*data->numN*sizeof(double));
    Jqp = (double*)malloc(data->numN*data->numN*sizeof(double));
    Jqq = (double*)malloc(data->numN*data->numN*sizeof(double));
    Pn = (double*)malloc(data->numN*sizeof(double));
    Qn = (double*)malloc(data->numN*sizeof(double));
    JacR = (double*)malloc(220*220*sizeof(double));

    while (Error<=100){
        calcularJacobiano(data,ybusReal,ybusImag,Vn,An,Jpp,Jpq,Jqp,Jqq,Pn,Qn);
        int i,j;
        for (i = 0 ; i < NumP ; i++) {
            N1 = NNP[i] - 1;
            dP[i] = Pref[N1] - Pn[N1];
        }

        for (i = 0; i < NumQ; i++ ) {
            N1 = NNQ[i] - 1;
            dQ[i] = Qref[N1] - Qn[N1];
        }

        createJacR(NNP, NNQ, NumQ, NumP, (int)data->numN, Jpp, Jpq, Jqp, Jqq, JacR);

        for (i = 0; i < NumP; i++) {
           for (j = 0; j < NumP; j++) {
               if(j!=NumP-1)
                   printf("%.4lf ",JacR[i*NumP+j]);
               else
                   printf("%.4lf\n",JacR[i*NumP+j]);

            }


            //printf("%.4lf ",dP[i]);
        }

        Error++;
    }

    /*int i;
    for (i = 0; i < data->numN; i++) {
        printf("%.4lf %.4lf\n",Pref[i],Qref[i]);
    }*/
   // printf("%.5lf\n",Vn[100]);
    //printf("%.5lf\n",data->numN);
    //printData(data,widthLineas, heightLineas, widthCargas, heightCargas, widthGen, heightGen);
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
    return res;
}
