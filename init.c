#include <stdio.h>
#include "utilities/utilities.h"

int main(){
    int res;
    char *fileName = "lineas";
    structData *data;
    data = NULL;
    res = loadDataFromFile(fileName, data);
    return res;
}
