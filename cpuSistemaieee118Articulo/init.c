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
    int mi = 186, n = 118, lda = 186, incx = 1, incy = 1;
    double alpha = 1, beta = 0, *Vrama, *Irama, *VnReal, *VnImag;

    int ni = 100, m;
    int k;
    char nameVrama[10], nameIrama[10], nameZpReal[10];
    double r, lambda = 0.0, kk, N, Pmax, Vw, Pw;
    NumW = 3;
    double Vmin = 4, Vnom = 12, Vmax = 25;
    srand((unsigned)time(NULL));
    Vrama = malloc(186*sizeof(double));
    Irama = malloc(186*sizeof(double));
    VnReal = malloc(118*sizeof(double));
    VnImag = malloc(118*sizeof(double));
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
        res = newtonRaphson(data, Vn);
        calculoVn(Vn, heightLineas, data->numN, VnReal, VnImag);
        //sprintf(nameVrama, "Vrama%d",k);
        //sprintf(nameIrama, "Irama%d",k);
        //sprintf(nameZpReal, "Zp%d",k);
       // printDataToFileVec(name,data->numN,Vn);
        dgemv_("N",&mi, &n, &alpha, A, &lda, Vn, &incx, &beta, Vrama, &incy);
        //printDataToFileVec(nameVrama,186,Vrama);
        //calculoIrama(Vrama,ZpReal,ZpImag,heightLineas,Irama);
        //printDataToFileVec(nameIrama,186,Irama);
        //printDataToFileVec(nameZpReal,186,ZpReal);




    }
    /////////////////////////////////////////////////////////////////////////////////////

    free(Imax);
    free(Irama);
    free(Vrama);
    free(data);
    free(Vn);
    return res;
}
