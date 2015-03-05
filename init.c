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
    structData *data;
    data = (structData*)malloc(sizeof(structData));
    res = loadDataFromFile(fileNameLineas,fileNameCargas,fileNameGen, data);
    printData(data,widthLineas, heightLineas, widthCargas, heightCargas, widthGen, heightGen);
    free(data);
    return res;
}
