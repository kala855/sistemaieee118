#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include "utilities/utilities.h"
#include <cusolverDn.h>
#include "utilities/cusolverUtilities.cuh"

/********************/
/* CUDA ERROR CHECK */
/********************/
// --- Credit to http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) { exit(code); }
    }
}

void gpuErrchk(cudaError_t ans) { gpuAssert((ans), __FILE__, __LINE__); }



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
    double *d_Jpp, *d_Jpq, *d_Jqp, *d_Jqq, *d_Pn, *d_Qn, *d_ybusReal, *d_ybusImag, *d_Vn, *d_An;
    double *d_dX, *d_JacRt, *d_work, *d_dP, *d_dQ, *d_Pref, *d_Qref, *d_JacR;
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

    int blockSize = 32;
    dim3 dimBlock(blockSize,1,1);
    dim3 dimGrid(ceil(data->numN/float(blockSize)),1,1);

    dim3 dimBlock2(blockSize,blockSize,1);
    dim3 dimGrid2(ceil(data->numN/float(blockSize)),ceil(data->numN/float(blockSize)),1);

    dim3 dimGrid3(ceil((data->numN*data->numN)/float(blockSize)),1,1);
    dim3 dimGrid4(ceil(NumP/float(blockSize)),1,1);
    dim3 dimGrid5(ceil(NumQ/float(blockSize)),1,1);
    dim3 dimGrid6(ceil(NumP/float(blockSize)), ceil(NumP/float(blockSize)));
    dim3 dimGrid7(ceil(NumP/float(blockSize)),1,1);

    t = 0;
    while (Error>1e-8){

        /*---- Initialize d_Jpp, d_Jpq, d_Jqp, d_Jqq, ----*/
        d_zeros<<<dimGrid3,dimBlock>>>(data->numN*data->numN,d_Jpp);
        d_zeros<<<dimGrid3,dimBlock>>>(data->numN*data->numN,d_Jpq);
        d_zeros<<<dimGrid3,dimBlock>>>(data->numN*data->numN,d_Jqp);
        d_zeros<<<dimGrid3,dimBlock>>>(data->numN*data->numN,d_Jqq);
        cudaDeviceSynchronize();

        gpuErrchk(cudaMemcpy(d_Vn,Vn,sizeof(double)*data->numN,cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(d_An,An,sizeof(double)*data->numN,cudaMemcpyHostToDevice));

        d_calcularJacobiano_1<<<dimGrid,dimBlock>>>(data->numN, d_ybusReal, d_ybusImag,d_Vn,\
                d_An,d_Pn,d_Qn);
        cudaDeviceSynchronize();
        d_calcularJacobiano_2<<<dimGrid2,dimBlock2>>>(data->numN, d_ybusReal, d_ybusImag, d_Vn\
                ,d_An, d_Pn,d_Qn, d_Jpp, d_Jpq, d_Jqp, d_Jqq);
        cudaDeviceSynchronize();

        gpuErrchk(cudaMemcpy(Pn,d_Pn,sizeof(double)*data->numN,cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(Qn,d_Qn,sizeof(double)*data->numN,cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(Jpp,d_Jpp,sizeof(double)*data->numN*data->numN,cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(Jpq,d_Jpq,sizeof(double)*data->numN*data->numN,cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(Jqp,d_Jqp,sizeof(double)*data->numN*data->numN,cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(Jqq,d_Jqq,sizeof(double)*data->numN*data->numN,cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(An,d_An,sizeof(double)*data->numN,cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(Vn,d_Vn,sizeof(double)*data->numN,cudaMemcpyDeviceToHost));

        if(t==0)
        {
            printDataToFileVec("pnData",data->numN,Pn);
            printDataToFileMat("jppData",data->numN,Jpp);
            printDataToFileMat("jpqData",data->numN,Jpq);
            printDataToFileMat("jqpData",data->numN,Jqp);
            printDataToFileMat("jqqData",data->numN,Jqq);
            printDataToFileVec("qnData",data->numN,Qn);
            t=1;
        }

        dp_compute<<<dimGrid4,dimBlock>>>(NumP, d_NNP, d_Pref, d_Pn, d_dP);
        gpuErrchk(cudaMemcpy(dP,d_dP,sizeof(double)*NumP,cudaMemcpyDeviceToHost));

        dq_compute<<<dimGrid5,dimBlock>>>(NumP, d_NNQ, d_Qref, d_Qn, d_dQ);
        gpuErrchk(cudaMemcpy(dQ,d_dQ,sizeof(double)*NumQ,cudaMemcpyDeviceToHost));

        //Se esta copiando JacR del Host al device. OJO
        gpuErrchk(cudaMemcpy(d_JacR,JacR,sizeof(double)*NumPQ*NumPQ,cudaMemcpyHostToDevice));
        d_createJacR_1<<<dimGrid6,dimBlock2>>>(d_NNP, NumQ, NumP, (int)(data->numN), d_Jpp, d_JacR);
        cudaDeviceSynchronize();
        d_createJacR_2<<<dimGrid7,dimBlock>>>(d_NNP, d_NNQ, NumQ, NumP, (int)(data->numN), d_Jpq, \
                d_JacR);
        cudaDeviceSynchronize();
        d_createJacR_3<<<dimGrid7,dimBlock>>>(d_NNP, d_NNQ, NumQ, NumP,(int)(data->numN), d_Jqp, \
                d_JacR);
        cudaDeviceSynchronize();
       /*d_createJacR_4<<<dimGrid5,dimBlock>>>(d_NNQ, NumQ, NumP,(int)(data->numN),\
                d_Jqq, d_JacR);
        cudaDeviceSynchronize();*/
        gpuErrchk(cudaMemcpy(JacR,d_JacR,sizeof(double)*NumPQ*NumPQ,cudaMemcpyDeviceToHost));


        createJacR(NNP, NNQ, NumQ, NumP, (int)data->numN, Jpp, Jpq, Jqp, Jqq, JacR);
        transposeJacR(JacR,NumPQ,JacRt);
        createdPdQ(dP,dQ,NumP,NumQ,dPdQ);
        cudaMemcpy(d_dX,dPdQ,sizeof(double)*NumPQ,cudaMemcpyHostToDevice);
        cudaMemcpy(d_JacRt,JacRt,sizeof(double)*NumPQ*NumPQ,cudaMemcpyHostToDevice);

        cusolveSafeCall(cusolverDnDgetrf(handle, NumPQ, NumPQ,d_JacRt,NumPQ,d_work,devIpiv,devInfo));
        gpuErrchk(cudaDeviceSynchronize());
        cusolveSafeCall(cusolverDnDgetrs(handle,trans,NumPQ,nrhs,d_JacRt,NumPQ,devIpiv,d_dX,NumPQ,devInfo));
        gpuErrchk(cudaMemcpy(dX,d_dX,sizeof(double)*NumPQ,cudaMemcpyDeviceToHost));



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
    // printDataToFileVec("pnData",data->numN,Pn);
    // printDataToFileVec("qnData",data->numN,Qn);
    printDataToFileMat("ybusRealData",data->numN,ybusReal);
    printDataToFileMat("ybusImagData",data->numN,ybusImag);

       /* printDataToFileMat("jppData",data->numN,Jpp);
          printDataToFileMat("jpqData",data->numN,Jpq);
          printDataToFileMat("jqpData",data->numN,Jqp);
          printDataToFileMat("jqqData",data->numN,Jqq);*/


    cusolverDnDestroy(handle);
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
    return res;
}
