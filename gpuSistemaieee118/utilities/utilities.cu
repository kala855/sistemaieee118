#include <stdlib.h>
#include "utilities.h"
#include <string.h>
#include <math.h>


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



__global__ void d_calc_An(double *dX, int *NNP, int NumP, double *An){
    int k = blockIdx.x*blockDim.x+threadIdx.x;
    int N1;
    if (k<NumP) {
        N1 = NNP[k] - 1;
        An[N1] = An[N1] + dX[k];
    }
}

__global__ void d_calc_Vn(double *dX, int *NNQ, int NumP, int NumQ, double *Vn){
    int k = blockIdx.x*blockDim.x+threadIdx.x;
    int kk, N1;
    if (k<NumQ) {
        N1 = NNQ[k] - 1;
        kk = k + NumP;
        Vn[N1] = Vn[N1] + dX[kk];
    }
}


__global__ void d_fill_d_dx(double *d_dPdQ, int NumPQ, double *d_dX){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<NumPQ){
        d_dX[i] = d_dPdQ[i];
    }
}


__global__ void d_filldPdQ1(double *d_dP, int NumP, double *d_dPdQ){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<NumP){
        d_dPdQ[i] = d_dP[i];
    }
}



__global__ void d_filldPdQ(double *d_dQ, int NumQ, int NumP, double *d_dPdQ){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<NumQ){
        d_dPdQ[i+NumP] = d_dQ[i];
    }
}

__global__ void d_transposeJacr(double *JacR,int NumPQ, double *JacRt){
    int i = blockIdx.y*blockDim.y+threadIdx.y;
    int j = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<NumPQ && j <NumPQ)
        JacRt[i*NumPQ+j] = JacR[j*NumPQ+i];

}

__global__ void d_zeros2(int size,double *An){
    int i = blockIdx.y*blockDim.y+threadIdx.y;
    int j = blockIdx.x*blockDim.x+threadIdx.x;
    if(i < size && j < size) {
        An[i*size+j] = 0.0;
    }
}


__global__ void d_zeros(int size,double *An){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i < size) {
        An[i] = 0.0;
    }
}

__global__ void d_ones(int size,double *Vn){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if( i < size){
        Vn[i] = 1.0;
    }
}

__global__ void dp_compute(int NumP,int *NNP, double *Pref, double *Pn, double *dP){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int N1;
    if(i<NumP){
        N1 = NNP[i] - 1;
        dP[i] = Pref[N1] - Pn[N1];
    }
}

__global__ void dq_compute(int NumQ,int *NNQ, double *Qref, double *Qn,double *dQ){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int N1;
    if(i<NumQ){
        N1 = NNQ[i] - 1;
        dQ[i] = Qref[N1] - Qn[N1];
    }
}



__global__ void d_createJacR_1(int *NNP, int NumQ, int NumP, int numN, double *Jpp, double *JacR){
    int i = blockIdx.y*blockDim.y+threadIdx.y;
    int j = blockIdx.x*blockDim.x+threadIdx.x;
    int size = NumP+NumQ;
    if((i<NumP) && (j<NumP)){
        JacR[i*size+j] =  Jpp[(NNP[i]-1)*numN+(NNP[j]-1)];
    }
}


__global__ void d_createJacR_2(int *NNP, int *NNQ, int NumQ, int NumP, int numN, double *Jpq, \
        double *JacR){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int size = NumP+NumQ;
    int j = 0;
    int k = i,l;
    if(i<NumP){
        l = NumP;
        for(j=0; j<NumQ; j++){
            JacR[k*size+l] = Jpq[(NNP[i]-1)*numN+(NNQ[j]-1)];
            l++;
        }
        //k++;
    }
}



__global__ void d_createJacR_3(int *NNP, int *NNQ, int NumQ, int NumP, int numN,double *Jqp, \
        double *JacR){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = 0;
    int k = NumP+i,l;
    int size = NumP+NumQ;
    if(i<NumQ){
        l = 0;
        for(j=0; j<NumP; j++){
            JacR[k*size+l] = Jqp[(NNQ[i]-1)*numN+(NNP[j]-1)];
            l++;
        }
        //k++;
    }
}



__global__ void d_createJacR_4(int *NNQ, int NumQ, int NumP, int numN, double *Jqq, double *JacR){
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int size = NumP+NumQ;
    int j = 0;
    int k = NumP+i,l;
    if(i<NumQ){
        l = NumP;
        for(j=0; j<NumQ; j++){
            JacR[k*size+l] = Jqq[(NNQ[i]-1)*numN+(NNQ[j]-1)];
            l++;
        }
    }
}

/*Esta funcion permite actualizar los valores de Pn y Qn que posteriormente serán usados
 para el cálculo completo del Jacobiano*/
__global__ void d_calcularJacobiano_1(int numN, double *ybusReal, double *ybusImag, double *Vn, \
        double *An, double *Pn, double *Qn){
    int k = blockIdx.x*blockDim.x+threadIdx.x;
    int widthLineas = numN;
    int m;
    double akm;
    if (k<numN){
        Pn[k] = 0.0;
        Qn[k] = 0.0;
        for (m = 0; m < numN; m++){
            akm = An[k] - An[m];
            Pn[k] = Pn[k] + Vn[m]*(ybusReal[k*widthLineas+m]*cos(akm) + \
                    ybusImag[k*widthLineas+m]*sin(akm));
            Qn[k] = Qn[k] - Vn[m]*(ybusImag[k*widthLineas+m]*cos(akm) - \
                    ybusReal[k*widthLineas+m]*sin(akm));

        }
        Pn[k] = Pn[k]*Vn[k];
        Qn[k] = Qn[k]*Vn[k];
    }


}

__global__ void d_calcularJacobiano_2(int numN, double *ybusReal, double *ybusImag, double *Vn,\
        double *An,double *Pn, double *Qn, double *Jpp, double *Jpq, double *Jqp, double *Jqq){

    int k = blockIdx.y*blockDim.y+threadIdx.y;
    int m = blockIdx.x*blockDim.x+threadIdx.x;
    int widthLineas = numN;
    double akm;
    if(k<numN){
        if(m<numN)
            /*Jpp[k*widthLineas+m] = 0.0;
            Jpq[k*widthLineas+m] = 0.0;
            Jqp[k*widthLineas+m] = 0.0;
            Jqq[k*widthLineas+m] = 0.0;*/

            if(k==m){
                Jpp[k*widthLineas+k] = -ybusImag[k*widthLineas+k]*Vn[k]*Vn[k]-Qn[k];
                Jpq[k*widthLineas+k] = ybusReal[k*widthLineas+k]*Vn[k] + Pn[k]/Vn[k];
                Jqp[k*widthLineas+k] = -ybusReal[k*widthLineas+k]*Vn[k]*Vn[k] + Pn[k];
                Jqq[k*widthLineas+k] = -ybusImag[k*widthLineas+k]*Vn[k] + Qn[k]/Vn[k];
            }else{
                akm = An[k] - An[m];
                Jpp[k*widthLineas+m] = Vn[k]*Vn[m]*(ybusReal[k*widthLineas+m]*sin(akm)- \
                        ybusImag[k*widthLineas+m]*cos(akm)) ;
                Jpq[k*widthLineas+m] = Vn[k]*(ybusReal[k*widthLineas+m]*cos(akm)+ybusImag[k*widthLineas+m]*sin(akm));
                Jqp[k*widthLineas+m] = -Jpq[k*widthLineas+m]*Vn[m];
                Jqq[k*widthLineas+m] = Jpp[k*widthLineas+m]/Vn[m];
            }
    }

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

int loadDataFromFile(char *filenameLineas, char *filenameCargas, char *filenameGen, structData *data){
    FILE *datosLineas,*datosGen,*datosCargas;
    int i, j, widthLineas,heightLineas,widthGen, heightGen,widthCargas, heightCargas;
    widthLineas = 6;
    heightLineas = 1896;/*186*/
    widthGen = 3;
    heightGen = 15*10;
    widthCargas = 3;
    heightCargas = 83*10;

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

    data->lineas = (double*)malloc(widthLineas*heightLineas*sizeof(double));
    data->cargas = (double*)malloc(widthCargas*heightCargas*sizeof(double));
    data->gen = (double*)malloc(widthGen*heightGen*sizeof(double));

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
