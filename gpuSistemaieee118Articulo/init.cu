#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utilities/utilities.h"
#include <time.h>
#include <math.h>
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,\
        double* b, int* ldb, int* info );

extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,\
        int *lda, double *x, int *incx, double *beta, double *y, int *incy);


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
    double *Imax, *A, *ZpReal, *ZpImag, *NW, *Vn, *An;
    structData *data;
    Mont *mont;
    Imax = malloc(numDataImax*sizeof(double));
    data = (structData*)malloc(sizeof(structData));
    mont = (Mont*)malloc(sizeof(Mont));
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
    Vn = (double*)malloc(data->numN*sizeof(double));
    An = (double*)malloc(data->numN*sizeof(double));
    ones(data->numN,Vn);
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


    //dgemv_ configuration parameters
    int mi = 186, n = 118, lda = 186, ldan = 118, incx = 1, incy = 1;
    double alpha = 1, beta = 0, *Vrama, *Irama, *VnReal, *VnImag, *ybusReal;
    double *ybusImag, *InReal, *InImag, *InImagTmp, *IlineaReal, *IlineaImag, *sobrecarga, *Vmedia;
    double *Vdesv, *Probmin, *Probmax, *Probsobrecarga;

    int ni = 1000, m;
    int k;
    char nameVrama[10], nameIrama[10], nameZpReal[10], nameVnReal[10], nameVnImag[10];
    char nameInReal[10], nameInImag[10], nameSobrecarga[10];
    double r, lambda = 0.0, kk, N, Pmax, Vw, Pw;
    NumW = 3;
    double Vmin = 4, Vnom = 12, Vmax = 25;
    srand((unsigned)time(NULL));
    ybusReal = malloc(data->numN*data->numN*sizeof(double));
    ybusImag = malloc(data->numN*data->numN*sizeof(double));
    InReal = malloc(data->numN*sizeof(double));
    InImag = malloc(data->numN*sizeof(double));
    InImagTmp = malloc(data->numN*sizeof(double));
    Vrama = malloc(data->numL*sizeof(double));
    Irama = malloc(data->numL*sizeof(double));
    IlineaReal = malloc(data->numL*sizeof(double));
    IlineaImag = malloc(data->numL*sizeof(double));
    sobrecarga = malloc(data->numL*sizeof(double));
    VnReal = malloc(data->numN*sizeof(double));
    VnImag = malloc(data->numN*sizeof(double));
    mont->sum = malloc(data->numN*sizeof(double));
    mont->sumcuad = malloc(data->numN*sizeof(double));
    mont->lv = malloc(data->numN*sizeof(double));
    mont->hv = malloc(data->numN*sizeof(double));
    mont->sob = malloc(data->numL*sizeof(double));
    Vmedia = malloc(data->numN*sizeof(double));
    Vdesv = malloc(data->numN*sizeof(double));
    Probmin = malloc(data->numN*sizeof(double));
    Probmax = malloc(data->numN*sizeof(double));
    Probsobrecarga = malloc(data->numL*sizeof(double));

    //mont Inicializacion
    zeros(data->numN, mont->sum);
    zeros(data->numN, mont->sumcuad);
    zeros(data->numN, mont->lv);
    zeros(data->numN, mont->hv);
    zeros(data->numL, mont->sob);

    for (k = 0; k < ni; k++) {
        for (m = 0; m < NumW; m++) {
            r = ((double)rand()/(double)RAND_MAX);
            lambda = NW[m*NumW+2];
            kk = NW[m*NumW+3];
            Pmax = NW[m*NumW+1];
            N = NW[m*NumW+0] - 1.0;
            Vw = lambda*pow((-log(1-r)),(1/kk));
            if ((Vw<Vmin)||(Vw>Vmax))
                Pw = 0.0;
            if ((Vmin<Vw)&&(Vw<Vnom))
                Pw = Pmax*pow((Vw/Vnom),3);
            if ((Vw>Vnom)&&(Vw<Vmax))
                Pw = Pmax;

            data->gen[(int)(N)*widthGen+1] = Pw;
            //printf("%lf \n",data->gen[(int)(N)*widthGen+1]);
        }
        // Llamar versión paralela
        res = newtonRaphson(data, Vn, An, ybusReal, ybusImag);
        calculoVn(Vn,An, data->numN, VnReal, VnImag);
        //cálculo de InReal --> Cálculo a Realizar en GPU
        dgemv_("N",&n, &n, &alpha, ybusReal, &ldan, VnReal, &incx, &beta, InReal, &incy);
        dgemv_("N",&n, &n, &alpha, ybusImag, &ldan, VnImag, &incx, &beta, InImag, &incy);
        res = subVectors(InReal, InImag, data->numN);

        //cálculo de InImag --> Cálculo a Realizar en GPU
        dgemv_("N",&n, &n, &alpha, ybusReal, &ldan, VnImag, &incx, &beta, InImag, &incy);
        dgemv_("N",&n, &n, &alpha, ybusImag, &ldan, VnReal, &incx, &beta, InImagTmp, &incy);
        res = addVectors(InImag,InImagTmp, data->numN);

        //cálculo de ilinea --> Cálculo a Realizar en GPU
        dgemv_("N",&mi, &n, &alpha, A, &lda, InReal, &incx, &beta, IlineaReal, &incy);
        dgemv_("N",&mi, &n, &alpha, A, &lda, InImag, &incx, &beta, IlineaImag, &incy);

        // Cálculo a Realizar en GPU
        calculoSobrecarga(IlineaReal, IlineaImag, sobrecarga, Imax, data->numL);
        calculoMontSobrecarga(data->numL, sobrecarga, mont);

        //Cálculo a Realizar en GPU
        calculoCorrientesRama(data->numN,Vn,mont);

  //      sprintf(nameSobrecarga, "sobrecarga%d",k);
        //sprintf(nameVnReal, "VnReal%d",k);
        //sprintf(nameVnImag, "VnImag%d",k);
        //sprintf(nameInReal, "InReal%d",k);
        //sprintf(nameInImag, "InImag%d",k);
        //sprintf(nameVrama, "Vrama%d",k);
        //sprintf(nameIrama, "Irama%d",k);
        //sprintf(nameZpReal, "Zp%d",k);
//        printDataToFileVec(nameSobrecarga,data->numL,mont->sob);
       // printDataToFileVec(name,data->numN,Vn);
        //printDataToFileVec(nameVnReal,data->numN,VnReal);
        //printDataToFileVec(nameVnReal,data->numN,VnReal);
        //printDataToFileVec(nameInReal,data->numL, IlineaReal);
        //printDataToFileVec(nameInImag,data->numL, IlineaImag);
        //dgemv_("N",&mi, &n, &alpha, A, &lda, Vn, &incx, &beta, Vrama, &incy);
        //printDataToFileVec(nameVrama,186,Vrama);
        //calculoIrama(Vrama,ZpReal,ZpImag,heightLineas,Irama);
        //printDataToFileVec(nameIrama,186,Irama);
        //printDataToFileVec(nameZpReal,186,ZpReal);
    }
    /////////////////////////////////////////////////////////////////////////////////////

    calculosFinales(data->numN, ni, mont,Vmedia,Vdesv,Probmin,Probmax);
    calculoProbSobrecarga(data->numL,ni,mont, Probsobrecarga);
    printDataToFileVec("vmedia",data->numN,Vmedia);
    printDataToFileVec("vdesv", data->numN,Vdesv);
    printDataToFileVec("probmin", data->numN,Probmin);
    printDataToFileVec("probmax", data->numN,Probmax);
    printDataToFileVec("probsobrecarga", data->numL,Probsobrecarga);

    free(Imax);
    free(VnReal);free(VnImag);
    free(ybusReal);free(ybusImag);
    free(InReal);free(InImag);free(InImagTmp);
    free(IlineaReal);free(IlineaImag);
    free(mont->sob);free(mont->lv);free(mont->sum);free(mont->sumcuad);free(mont->hv);free(mont);
    free(Vmedia);free(Vdesv);free(Probmin);free(Probmax);free(Probsobrecarga);
    free(sobrecarga);
    free(Irama);
    free(Vrama);
    free(data);
    free(Vn);
    return res;
}
