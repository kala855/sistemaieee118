#include <stdlib.h>
#include "utilities.h"
#include <string.h>

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
                                               - tap*tap*YYImag;
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
