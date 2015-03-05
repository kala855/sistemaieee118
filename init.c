#include <stdio.h>
#include <stdlib.h>
#include "utilities/utilities.h"

int main(){
    int res;
    int widthCargas = 3;
    int widthGen = 3;
    int widthLineas = 5;
    int heightCargas = 83;
    int heightGen = 15;
    int heightLineas = 186;
    char *fileNameLineas = "lineas";
    char *fileNameCargas = "cargas";
    char *fileNameGen = "gen";
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
    ones(data->numN,Vn);
    zeros(data->numN,An);
    printf("%.5lf",An[105]);
   // printf("%.5lf\n",Vn[100]);
    //printf("%.5lf\n",data->numN);
    //printData(data,widthLineas, heightLineas, widthCargas, heightCargas, widthGen, heightGen);
    free(data);
    free(Vn);
    free(An);
    return res;
}
