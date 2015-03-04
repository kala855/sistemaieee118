#include <stdlib.h>
#include "utilities.h"

int loadDataFromFile(char *filename, structData *data){
    FILE *datos;
    datos = fopen(filename,"r");
    if (datos == NULL){
        printf("Archivo de Lineas inexistente %s verifique \n", filename);
        exit(1);
    }
    data->lineas = malloc(5*186*sizeof(double));
    if (data->lineas == NULL){
        printf("Imposible asignar memoria a lineas\n");
        exit(1);
    }

    printf("test\n");
    //fscanf(datos, "%lf,%lf,%lf,%lf,%lf",&data->lineas[0],&data->lineas[1],&data->lineas[2],&data->lineas[3],&data->lineas[4]);

    //printf("%.5lf,%.5lf,%.5lf,%.5lf,%.5lf\n",data->lineas[0],data->lineas[1],data->lineas[2],data->lineas[3],data->lineas[4]);
    return 0;
    fclose(datos);
}
