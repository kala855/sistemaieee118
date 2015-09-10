#include <stdlib.h>
#include "utilities.h"
#include <cusolverUtilities.cuh>
#include <string.h>
#include <math.h>
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,\
        double* b, int* ldb, int* info );



int calcCargLineas(structData *data, double *An, double *Vn, double *Ism){
    int k, N1, N2;
    double akm;
    int widthLineas = 6;
    for (k = 0; k < data->numL; k++) {
        N1 = data->lineas[k*widthLineas+0] - 1;
        N2 = data->lineas[k*widthLineas+1] - 1;
        akm = An[N1] - An[N2];
        Ism[k] = sqrt(Vn[N1]*Vn[N1] + Vn[N2]*Vn[N2] - 2*Vn[N1]*Vn[N2]*cos(akm));
        Ism[k] = Ism[k]/sqrt(data->lineas[k*widthLineas+2]*data->lineas[k*widthLineas+2] + \
                data->lineas[k*widthLineas+3]*data->lineas[k*widthLineas+3])\
                 /1.6;//data->lineas[k*widthLineas+6];

    }
    return 0;
}

double maxAbs(int NumPQ, double *dPdQ){
    int i;
    double max = 0.0;
    for (i = 0; i < NumPQ; i++){
        if(fabs(dPdQ[i]) > max)
            max = fabs(dPdQ[i]);
    }
    return max;
}

int transposeJacR(double *JacR,int NumPQ, double *JacRt){
    int i, j;
    for ( i = 0; i < NumPQ; i++) {
        for (j = 0; j < NumPQ; j++) {
            JacRt[i*NumPQ+j] = JacR[j*NumPQ+i];
        }
    }
    return 0;

}

int createdPdQ(double *dP, double *dQ, int NumP, int NumQ, double *dPdQ){
    memcpy(dPdQ,dP,sizeof(double)*NumP);
    memcpy(dPdQ+NumP,dQ,sizeof(double)*NumQ);
    return 0;
}

int createJacR(int *NNP, int *NNQ, int NumQ, int NumP, int numN,double *Jpp, \
        double *Jpq, double *Jqp, double *Jqq, double *JacR){
   // printf("%d %d\n", NumQ, NumP);
    int size = NumP+NumQ;
    int i,j,k,l;
    for (i = 0; i < NumP; i++) {
        for (j = 0; j < NumP; j++) {
            JacR[i*size+j] = Jpp[(NNP[i]-1)*numN+(NNP[j]-1)];
        }
    }

    k = 0;
    for (i = 0; i < NumP; i++) {
        l = NumP;
        for (j = 0; j < NumQ; j++) {
            JacR[k*size+l] = Jpq[(NNP[i]-1)*numN+(NNQ[j]-1)];
            l++;
        }
        k++;
    }

    k = NumP;
    for (i = 0; i < NumQ; i++) {
        l = 0;
        for (j = 0; j < NumP; j++) {
            JacR[k*size+l] = Jqp[(NNQ[i]-1)*numN+(NNP[j]-1)];
            l++;
        }
        k++;
    }

    k = NumP;
    for (i = 0; i < NumQ; i++) {
        l = NumP;
        for (j = 0; j < NumQ; j++) {
            JacR[k*size+l] = Jqq[(NNQ[i]-1)*numN+(NNQ[j]-1)];
            l++;
        }
        k++;
    }



    return 0;
}

int calcularJacobiano(structData *data, double *ybusReal, double *ybusImag, double *Vn, double *An, \
        double *Jpp, double *Jpq, double *Jqp, double *Jqq, double *Pn, double *Qn){

    zeros(data->numN,Pn);
    zeros(data->numN,Qn);
    double akm, *H,*N,*J,*L;
    int k,m,widthLineas=data->numN;
    for (k = 0; k < data->numN; k++) {
        for (m = 0; m < data->numN; m++) {
            akm = An[k] - An[m];
            Pn[k] = Pn[k] + Vn[m]*(ybusReal[k*widthLineas+m]*cos(akm) + \
                    ybusImag[k*widthLineas+m]*sin(akm));
            Qn[k] = Qn[k] - Vn[m]*(ybusImag[k*widthLineas+m]*cos(akm) - \
                    ybusReal[k*widthLineas+m]*sin(akm));
        }
        Pn[k] = Pn[k]*Vn[k];
        Qn[k] = Qn[k]*Vn[k];
    }

    H = (double*)malloc(data->numN*data->numN*sizeof(double));
    N = (double*)malloc(data->numN*data->numN*sizeof(double));
    J = (double*)malloc(data->numN*data->numN*sizeof(double));
    L = (double*)malloc(data->numN*data->numN*sizeof(double));

    zeros(data->numN*data->numN,H);
    zeros(data->numN*data->numN,N);
    zeros(data->numN*data->numN,J);
    zeros(data->numN*data->numN,L);

    for (k = 0; k < data->numN; k++) {
        for (m = 0; m < data->numN; m++) {
            if(k==m){
                H[k*widthLineas+k] = -ybusImag[k*widthLineas+k]*Vn[k]*Vn[k]-Qn[k];
                N[k*widthLineas+k] = ybusReal[k*widthLineas+k]*Vn[k] + Pn[k]/Vn[k];
                J[k*widthLineas+k] = -ybusReal[k*widthLineas+k]*Vn[k]*Vn[k] + Pn[k];
                L[k*widthLineas+k] = -ybusImag[k*widthLineas+k]*Vn[k] + Qn[k]/Vn[k];
            }else{
                akm = An[k] - An[m];
                H[k*widthLineas+m] = Vn[k]*Vn[m]*(ybusReal[k*widthLineas+m]*sin(akm)- \
                        ybusImag[k*widthLineas+m]*cos(akm)) ;
                N[k*widthLineas+m] = Vn[k]*(ybusReal[k*widthLineas+m]*cos(akm)+ybusImag[k*widthLineas+m]*sin(akm));
                J[k*widthLineas+m] = -N[k*widthLineas+m]*Vn[m];
                L[k*widthLineas+m] = H[k*widthLineas+m]/Vn[m];
            }
        }
    }

    memcpy(Jpp,H,data->numN*data->numN*sizeof(double));
    memcpy(Jpq,N,data->numN*data->numN*sizeof(double));
    memcpy(Jqp,J,data->numN*data->numN*sizeof(double));
    memcpy(Jqq,L,data->numN*data->numN*sizeof(double));

    free(H);free(N);free(J);free(L);
    return 0;

}

int setdiff(int *vector1, double *vector2, int size1, int size2, int *c){
    int i, j, accum=0, k=0;
    for (i = 0; i < size1; i++) {
        for (j = 0; j < size2; j++) {
            if(vector1[i]!=(int)vector2[j*3])
                accum++;
            else
                break;
        }
        if(accum == size2){
            c[k] = vector1[i];
            k++;
        }
        accum = 0;
    }
    return k;
}

int genVector(int *NNP, int initNumber,int finalNumber){
    int i = initNumber,j=0;
    for (i = initNumber; i <= finalNumber; i++) {
        NNP[j] = i;
        j++;
    }
    return 0;
}

int calcularYbus(structData *data, double *ybusReal, double *ybusImag){
    zeros(data->numN*data->numN, ybusReal);
    zeros(data->numN*data->numN, ybusImag);
    int widthLineas = 6;
    double N1, N2, ym;
    double YYReal, YYImag, tap;
    int k,t;
    for (k = 0; k < data->numL; k++) {
        N1 = data->lineas[k*widthLineas+0]-1;
        N2 = data->lineas[k*widthLineas+1]-1;
        ym = data->lineas[k*widthLineas+2]*data->lineas[k*widthLineas+2]+\
             data->lineas[k*widthLineas+3]*data->lineas[k*widthLineas+3];
        YYReal = data->lineas[k*widthLineas+2]/ym;
        YYImag = -data->lineas[k*widthLineas+3]/ym;

        tap = 1/data->lineas[k*widthLineas+5];

        ybusReal[(int) (N1*data->numN + N1)] = ybusReal[(int) (N1*data->numN + N1)]\
                                               + tap*tap*YYReal;
        ybusReal[(int) (N2*data->numN + N2)] = ybusReal[(int) (N2*data->numN + N2)] \
                                               + YYReal;
        ybusReal[(int) (N1*data->numN + N2)] = ybusReal[(int) (N1*data->numN + N2)] \
                                               - tap*YYReal;
        ybusReal[(int) (N2*data->numN + N1)] = ybusReal[(int) (N2*data->numN + N1)] \
                                               - tap*YYReal;

        ybusImag[(int) (N1*data->numN + N1)] = ybusImag[(int) (N1*data->numN + N1)] \
                                               + tap*tap*YYImag + data->lineas[k*widthLineas+4];
        ybusImag[(int) (N2*data->numN + N2)] = ybusImag[(int) (N2*data->numN + N2)] \
                                               + YYImag + data->lineas[k*widthLineas+4];
        ybusImag[(int) (N1*data->numN + N2)] = ybusImag[(int) (N1*data->numN + N2)] \
                                               - tap*YYImag;
        ybusImag[(int) (N2*data->numN + N1)] = ybusImag[(int) (N2*data->numN + N1)] \
                                               - tap*YYImag;
    }

    return 0;


}

int zeros(int size,double *An){
    int i;
    for (i = 0; i < size; i++) {
        An[i] = 0.0;
    }
    return 0;
}

int ones(int size,double *Vn){
    int i;
    for (i = 0; i < size; i++) {
        Vn[i] = 1.0;
    }
    return 0;
}

double maxLineas(structData *data, int widthLineas,int heightLineas){
    int i;
    double mayorLineas1 = 0, mayorLineas2;
    for (i = 0; i < heightLineas; i++) {
        if(mayorLineas1 < data->lineas[i*widthLineas]){
            mayorLineas1 = data->lineas[i*widthLineas];
        }
        if(mayorLineas2 < data->lineas[i*widthLineas+1]){
            mayorLineas2 = data->lineas[i*widthLineas+1];
        }
    }

    if(mayorLineas2>=mayorLineas1)
        return mayorLineas2;
    else
        return mayorLineas1;

}

int printData(structData *data, int widthLineas, int heightLineas, int widthCargas, int heightCargas, int widthGen, int heightGen){
    int i;
    printf("Datos Líneas\n");
    for (i = 0; i < heightLineas; i++) {
        printf("%.5lf,%.5lf,%.5lf,%.5lf,%.5lf\n",data->lineas[i*widthLineas],data->lineas[i*widthLineas+1],data->lineas[i*widthLineas+2],data->lineas[i*widthLineas+3],data->lineas[i*widthLineas+4]);
    }

    printf("Datos Cargas\n");
    for (i = 0; i < heightCargas; i++) {
        printf("%.5lf,%.5lf,%.5lf\n",data->cargas[i*widthCargas],data->cargas[i*widthCargas+1],data->cargas[i*widthCargas+2]);
    }

    printf("Datos Gen\n");
    for (i = 0; i < heightGen; i++) {
        printf("%.5lf,%.5lf,%.5lf\n",data->gen[i*widthGen],data->gen[i*widthGen+1],data->gen[i*widthGen+2]);
    }

    return 0;
}

int loadCorrientesMax(char *fileNameIMax, double *Imax){
    FILE *datosImax;
    int numDataImax, i;
    numDataImax = 186;
    datosImax = fopen(fileNameIMax,"r");
    if (datosImax == NULL){
        printf("Archivo de Imax inexistente %s verifique \n",fileNameIMax);
        exit(1);
    }

    for(i=0;i<numDataImax;i++){
        fscanf(datosImax,"%lf\n",&Imax[i]);
    }

    fclose(datosImax);
    return 0;
}



int loadDataFromFile(char *filenameLineas, char *filenameCargas, char *filenameGen, structData *data){
    FILE *datosLineas,*datosGen,*datosCargas;
    int i, j, widthLineas,heightLineas,widthGen, heightGen,widthCargas, heightCargas;
    widthLineas = 6;
    heightLineas = 186;
    widthGen = 3;
    heightGen = 15;
    widthCargas = 3;
    heightCargas = 83;

    datosGen = fopen(filenameGen,"r");
    datosLineas = fopen(filenameLineas,"r");
    datosCargas = fopen(filenameCargas,"r");

    if (datosGen == NULL){
        printf("Archivo de Gen inexistente %s verifique \n", filenameGen);
        exit(1);
    }

    if (datosLineas == NULL){
        printf("Archivo de Lineas inexistente %s verifique \n", filenameLineas);
        exit(1);
    }

    if(datosCargas == NULL){
        printf("Archivo de Cargas inexistente %s verifique \n", filenameCargas);
        exit(1);
    }

    data->lineas = malloc(widthLineas*heightLineas*sizeof(double));
    data->cargas = malloc(widthCargas*heightCargas*sizeof(double));
    data->gen = malloc(widthGen*heightGen*sizeof(double));

    if (data->lineas == NULL){
        printf("Imposible asignar memoria a lineas\n");
        exit(1);
    }
    if (data->gen == NULL){
        printf("Imposible asignar memoria a gen\n");
        exit(1);
    }

    if (data->cargas == NULL){
        printf("Imposible asignar memoria a cargas\n");
        exit(1);
    }
    for (i = 0; i < heightLineas; i++) {//Se adiciona el tap constante en 1.0 para este caso
        fscanf(datosLineas, "%lf,%lf,%lf,%lf,%lf\n",&data->lineas[i*widthLineas],\
                &data->lineas[i*widthLineas+1],&data->lineas[i*widthLineas+2],\
                &data->lineas[i*widthLineas+3],&data->lineas[i*widthLineas+4]);
        if(data->lineas[i*widthLineas+4]==0)
            data->lineas[i*widthLineas+4] = 1.0;

        data->lineas[i*widthLineas+5] = 1.0;
    }

    for (i = 0; i < heightCargas; i++) {
        fscanf(datosCargas, "%lf,%lf,%lf\n",&data->cargas[i*widthCargas],&data->cargas[i*widthCargas+1],&data->cargas[i*widthCargas+2]);
    }

    for (i = 0; i < heightGen; i++) {
        fscanf(datosGen, "%lf,%lf,%lf\n",&data->gen[i*widthGen],&data->gen[i*widthGen+1],&data->gen[i*widthGen+2]);
    }

    fclose(datosLineas);
    fclose(datosGen);
    fclose(datosCargas);
    return 0;
}

int printDataToFileVec(char *name, int size,double *data){
    FILE *dato;
    dato = fopen(name,"w");
    int i;
    for (i = 0; i < size; i++) {
        fprintf(dato,"%.4lf\n",data[i]);
    }
    fclose(dato);
    return 0;

}

int printDataToFileMat(char *name, int size,double *data){
    FILE *dato;
    dato = fopen(name,"w");
    int i,j;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            if(j!=size-1)
                fprintf(dato,"%.4lf ",data[i*size+j]);
            else
                fprintf(dato,"%.4lf\n",data[i*size+j]);
        }

    }

    fclose(dato);
    return 0;

}

int calcularMatrizA(structData *data, int widthLineas, double *A){
    zeros((data->numL)*(data->numN),A);
    int i;
    int N1,N2;
    for(i=0;i<(data->numL);i++){
        N1 = (int)(data->lineas[i*widthLineas+0])-1;
        N2 = (int)(data->lineas[i*widthLineas+1])-1;
        A[i*(int)(data->numN)+N1] = 1;
        A[i*(int)(data->numN)+N2] = -1;
    }
    return 0;
}

int printMatrixToFile(double *A, int numFilas, int numColumnas, char *name){
    FILE *dato;
    dato = fopen(name,"w");
    int i,j;
    for (i = 0; i < numFilas; i++) {
        for (j = 0; j < numColumnas; j++) {
            if(j!=numColumnas-1)
                fprintf(dato,"%.0lf ", A[i*(int)(numColumnas)+j]);
            else
                fprintf(dato,"%.0lf\n", A[i*(int)(numColumnas)+j]);

        }
    }
    fclose(dato);
    return 0;
}

int calcularZp(structData *data, int heightLineas, int widthLineas, double *ZpReal, double *ZpImag){
    int i;
    for (i = 0; i < heightLineas; i++) {
        ZpReal[i] = data->lineas[i*widthLineas+2];
        ZpImag[i] = data->lineas[i*widthLineas+3];

    }
    return 0;
}


int loadNW(char *fileNameNW, double *NW){
    FILE *datosNW;
    int numDataNW, i;
    numDataNW = 3;
    datosNW = fopen(fileNameNW,"r");
    if (datosNW == NULL){
        printf("Archivo de NW inexistente %s verifique \n",fileNameNW);
        exit(1);
    }

    for(i=0;i<numDataNW;i++){
        fscanf(datosNW,"%lf,%lf,%lf,%lf\n",&NW[i*4+0],&NW[i*4+1],&NW[i*4+2],&NW[i*4+3]);
    }

    fclose(datosNW);
}

int newtonRaphsonCUDA(structData *data, double *Vn, double *ybusReal, double *ybusImag){
    int res,i,j;
    int widthCargas = 3;
    int widthGen = 3;
    int widthLineas = 6;
    int heightCargas = 83;
    int heightGen = 15;
    int heightLineas = 186, NumP;
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
    double *d_Jpp, *d_Jpq, *d_Jqp, *d_Jqq, *d_Pn, *d_Qn, *d_ybusReal, *d_ybusImag, *d_Vn, *d_An;
    double *d_dX, *d_JacRt, *d_work, *d_dP, *d_dQ, *d_Pref, *d_Qref, *d_JacR, *d_dPdQ;
    int *devIpiv, *d_NNP, *d_NNQ;

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

    cudaError_t error = cudaSuccess;
    int lwork =0, *devInfo;
    gpuErrchk(cudaMalloc(&devInfo, sizeof(int)));


    gpuErrchk(cudaMalloc(&d_ybusImag,data->numN*data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_ybusReal,data->numN*data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Qn,data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Pn,data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Jpp,data->numN*data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Jqp,data->numN*data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Jpq,data->numN*data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Jqq,data->numN*data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Vn,data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_An,data->numN*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_dP,NumP*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_dQ,NumQ*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Pref,(data->numN)*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_Qref,(data->numN)*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_NNP,((data->numN)-1)*sizeof(int)));
    gpuErrchk(cudaMalloc(&d_NNQ,(data->numN)*sizeof(int)));
    gpuErrchk(cudaMalloc(&d_JacR,NumPQ*NumPQ*sizeof(double)));
    gpuErrchk(cudaMalloc(&d_dPdQ,NumPQ*sizeof(double)));


    // ---- Copy ybusData to GPU ----//
    gpuErrchk(cudaMemcpy(d_ybusReal,ybusReal,sizeof(double)*data->numN*data->numN\
                ,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_ybusImag,ybusImag,sizeof(double)*data->numN*data->numN\
                ,cudaMemcpyHostToDevice));
    /*---- Copy Pref and Qref to device ----*/
    gpuErrchk(cudaMemcpy(d_Pref,Pref,sizeof(double)*data->numN,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_Qref,Qref,sizeof(double)*data->numN,cudaMemcpyHostToDevice));

    /*---- Copy Pref and Qref to device ----*/
    gpuErrchk(cudaMemcpy(d_NNP,NNP,sizeof(int)*((data->numN)-1),cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_NNQ,NNQ,sizeof(int)*data->numN,cudaMemcpyHostToDevice));
    /*---- Copy An y Vn al device -----*/
    gpuErrchk(cudaMemcpy(d_Vn,Vn,sizeof(double)*data->numN,cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_An,An,sizeof(double)*data->numN,cudaMemcpyHostToDevice));


    // ---- cuSolver initialization ---- //
    cusolverStatus_t solvStatus = CUSOLVER_STATUS_SUCCESS;
    cusolverDnHandle_t handle;
    solvStatus = cusolverDnCreate(&handle);
    ///////////////////////////////////////

    cublasOperation_t trans = CUBLAS_OP_N;


    gpuErrchk(cudaMalloc((void**)&d_JacRt,sizeof(double)*NumPQ*NumPQ));
    gpuErrchk(cudaMalloc((void**)&d_dX,sizeof(double)*NumPQ));
    gpuErrchk(cudaMalloc((void**)&devIpiv,sizeof(int)*NumPQ));

    cusolveSafeCall(cusolverDnDgetrf_bufferSize(handle,NumPQ,NumPQ,d_JacRt,lda,&lwork));

    gpuErrchk(cudaMalloc((void**)&d_work, sizeof(double)*lwork));

    zeros(data->numL,Ism);

    int blockSize2D = 32;
    int blockSize = 1024;
    dim3 dimBlock(blockSize,1,1);
    dim3 dimGrid(ceil(data->numN/float(blockSize)),1,1);

    dim3 dimBlockXY(blockSize2D,blockSize2D,1);
    dim3 dimGridXY(ceil(data->numN/float(blockSize2D)),ceil(data->numN/float(blockSize2D)),1);

    dim3 dimBlock2(blockSize2D,blockSize2D,1);
    dim3 dimGrid2(ceil(data->numN/float(blockSize2D)),ceil(data->numN/float(blockSize2D)),1);

    dim3 dimGrid3(ceil((data->numN*data->numN)/float(blockSize)),1,1);
    dim3 dimGrid4(ceil(NumP/float(blockSize)),1,1);
    dim3 dimGrid5(ceil(NumQ/float(blockSize)),1,1);
    dim3 dimGrid6(ceil(NumP/float(blockSize2D)), ceil(NumP/float(blockSize2D)),1);
    dim3 dimGrid7(ceil(NumP/float(blockSize)),1,1);
    dim3 dimGrid8(ceil(NumPQ/float(blockSize2D)),ceil(NumPQ/float(blockSize2D)),1);
    dim3 dimGrid9(ceil(NumQ/float(blockSize)),1,1);
    dim3 dimGrid10(ceil(NumPQ/float(blockSize)),1,1);

    while (Error>1e-8){

        /*---- Initialize d_Jpp, d_Jpq, d_Jqp, d_Jqq, ----*/
        d_zeros<<<dimGrid3,dimBlock>>>(data->numN*data->numN,d_Jpp);
        d_zeros<<<dimGrid3,dimBlock>>>(data->numN*data->numN,d_Jpq);
        d_zeros<<<dimGrid3,dimBlock>>>(data->numN*data->numN,d_Jqp);
        d_zeros<<<dimGrid3,dimBlock>>>(data->numN*data->numN,d_Jqq);
        cudaDeviceSynchronize();

        d_calcularJacobiano_1<<<dimGrid,dimBlock>>>(data->numN, d_ybusReal, d_ybusImag,d_Vn,\
                d_An,d_Pn,d_Qn);

        cudaDeviceSynchronize();
        d_calcularJacobiano_2<<<dimGrid2,dimBlock2>>>(data->numN, d_ybusReal, d_ybusImag, d_Vn\
                ,d_An, d_Pn,d_Qn, d_Jpp, d_Jpq, d_Jqp, d_Jqq);
        cudaDeviceSynchronize();

        dp_compute<<<dimGrid4,dimBlock>>>(NumP, d_NNP, d_Pref, d_Pn, d_dP);

        dq_compute<<<dimGrid5,dimBlock>>>(NumP, d_NNQ, d_Qref, d_Qn, d_dQ);

        d_createJacR_1<<<dimGrid6,dimBlock2>>>(d_NNP, NumQ, NumP, (int)(data->numN), d_Jpp, d_JacR);
        cudaDeviceSynchronize();
        d_createJacR_2<<<dimGrid7,dimBlock>>>(d_NNP, d_NNQ, NumQ, NumP, (int)(data->numN), d_Jpq, \
                d_JacR);
        cudaDeviceSynchronize();
        d_createJacR_3<<<dimGrid7,dimBlock>>>(d_NNP, d_NNQ, NumQ, NumP,(int)(data->numN), d_Jqp, \
                d_JacR);
        cudaDeviceSynchronize();
        d_createJacR_4<<<dimGrid5,dimBlock>>>(d_NNQ, NumQ, NumP,(int)(data->numN),\
                d_Jqq, d_JacR);
        cudaDeviceSynchronize();


        d_transposeJacr<<<dimGrid8,dimBlock2>>>(d_JacR, NumPQ, d_JacRt);
        cudaDeviceSynchronize();

        d_filldPdQ1<<<dimGrid7,dimBlock>>>(d_dP,NumP,d_dPdQ);
        cudaDeviceSynchronize();
        d_filldPdQ<<<dimGrid9,dimBlock>>>(d_dQ,NumQ,NumP,d_dPdQ);
        cudaDeviceSynchronize();

        d_fill_d_dx<<<dimGrid10,dimBlock>>>(d_dPdQ,NumPQ,d_dX);
        cudaDeviceSynchronize();


        cusolveSafeCall(cusolverDnDgetrf(handle, NumPQ, NumPQ,d_JacRt,NumPQ,d_work,devIpiv,devInfo));
        gpuErrchk(cudaDeviceSynchronize());
        cusolveSafeCall(cusolverDnDgetrs(handle,trans,NumPQ,nrhs,d_JacRt,NumPQ,devIpiv,d_dX,NumPQ,\
                    devInfo));
        cudaDeviceSynchronize();


        d_calc_An<<<dimGrid7,dimBlock>>>(d_dX, d_NNP, NumP, d_An);
        cudaDeviceSynchronize();

        d_calc_Vn<<<dimGrid9,dimBlock>>>(d_dX, d_NNQ, NumP, NumQ, d_Vn);
        cudaDeviceSynchronize();

        gpuErrchk(cudaMemcpy(dPdQ, d_dPdQ,sizeof(double)*NumPQ,cudaMemcpyDeviceToHost));
        Error = maxAbs(NumPQ,dPdQ);

        if (iter>data->maxIter) {
            printf("..... No converge despues de %d iteraciones\nError = %lf\n", data->maxIter, Error);
            break;
        }
        iter++;
    }


    gpuErrchk(cudaMemcpy(An,d_An,sizeof(double)*data->numN,cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(Vn,d_Vn,sizeof(double)*data->numN,cudaMemcpyDeviceToHost));
    calcCargLineas(data,An,Vn,Ism);
    printDataToFileVec("ismData",data->numL,Ism);
    printDataToFileVec("vnData",data->numN,Vn);
    printDataToFileVec("anData",data->numN,An);

    cusolverDnDestroy(handle);
    free(An);
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
    cudaFree(d_ybusReal);
    cudaFree(d_ybusImag);
    cudaFree(d_work);
    cudaFree(d_JacRt);
    cudaFree(devIpiv);
    cudaFree(d_work);
    cudaFree(d_Jpp);
    cudaFree(d_Jqq);
    cudaFree(d_Jpq);
    cudaFree(d_Jqp);
    cudaFree(d_Vn);
    cudaFree(d_An);
    cudaFree(d_Pref);
    cudaFree(d_Qref);
    cudaFree(d_dP);
    cudaFree(d_dQ);
    cudaFree(d_JacR);
    cudaFree(d_JacRt);
    cudaFree(d_dPdQ);
    return res;
}

int newtonRaphson(structData *data, double *Vn, double *An, double *ybusReal, double *ybusImag){
    int res,i,j;
    int widthCargas = 3;
    int widthGen = 3;
    int widthLineas = 6;
    int heightCargas = 83;
    int heightGen = 15;
    int heightLineas = 186, NumP;
    //Vn = (double*)malloc(data->numN*sizeof(double));
    //An = (double*)malloc(data->numN*sizeof(double));
//    double *ybusReal = (double*) malloc(data->numN*data->numN*sizeof(double));
  //  double *ybusImag = (double*) malloc(data->numN*data->numN*sizeof(double));
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
    //printDataToFileVec("ismData",data->numL,Ism);
    //printDataToFileVec("vnData",data->numN,Vn);
    //printDataToFileVec("anData",data->numN,An);

    //free(ybusReal);
    //free(ybusImag);
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

int calculoIrama(double *Vrama, double *ZpReal, double *ZpImag, int heightLineas, double *Irama){
    int i;
    for (i = 0; i < heightLineas; i++) {
        Irama[i] = fabs(Vrama[i]/ZpReal[i]);
    }
    return 0;
}


int calculoVn(double *Vn, double *An, int height, double *VnReal, double *VnImag){
    int i;
    for (i = 0; i < height; i++) {
        VnReal[i] = Vn[i]*cos(An[i]);
        VnImag[i] = Vn[i]*sin(An[i]);
    }
    return 0;
}

int addVectors(double *vec1, double *vec2, int numN){
    int i;
    for (i = 0; i < numN; i++) {
        vec1[i] = vec1[i] + vec2[i];
    }
    return 0;
}

int subVectors(double *vec1, double *vec2, int numN){
    int i;
    for (i = 0; i < numN; i++) {
        vec1[i] = vec1[i] - vec2[i];
    }
    return 0;
}

int calculoSobrecarga(double *IlineaReal, double *IlineaImag, double *sobrecarga, double *Imax, int numL){
    int i;
    for (i = 0; i < numL; i++) {
        sobrecarga[i] = sqrt(IlineaReal[i]*IlineaReal[i] + IlineaImag[i]*IlineaImag[i])/Imax[i];
    }
    return 0;
}

int calculoMontSobrecarga(int numL, double *sobrecarga, Mont *mont){
    int i;
    for (i = 0; i < numL; i++) {
        if(sobrecarga[i]>1.0)
            mont->sob[i] = mont->sob[i] + 1;
    }
    return 0;
}


int calculoCorrientesRama(int numN,double *Vn,Mont *mont){
    int i;
    for (i = 0; i < numN; i++) {
        mont->sum[i] = mont->sum[i] + Vn[i];
        mont->sumcuad[i] = mont->sumcuad[i] + Vn[i] * Vn[i];
        if(Vn[i]<0.9)
            mont->lv[i] = mont->lv[i] + 1;
        if(Vn[i]>1.1)
            mont->hv[i] = mont->hv[i] + 1;
    }
    return 0;
}

int calculosFinales(int numN,int ni, Mont *mont, double *Vmedia, double *Vdesv, double *Probmin, double *ProbMax){
    int i;
    for (i = 0; i < numN; i++) {
        Vmedia[i] = mont->sum[i]/ni;
        Vdesv[i] = mont->sumcuad[i] - 2*mont->sum[i]*Vmedia[i] - Vmedia[i]*Vmedia[i];
        Probmin[i] = mont->lv[i]/ni;
        ProbMax[i] = mont->hv[i]/ni;
    }
}

int calculoProbSobrecarga(int numL, int ni, Mont *mont, double *Probsobrecarga){
    int i;
    for (i = 0; i < numL; i++) {
        Probsobrecarga[i] = mont->sob[i]/ni;
    }
    return 0;
}
