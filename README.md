# sistemaieee118

### Compilación
El proyecto hace uso de cmake. Para realizar el proceso de compilación se debe hacer lo siguiente:

    git clone https://github.com/kala855/sistemaieee118.git

Con este paso se obtendrán todos los archivos fuente y de entrada que necesita el proyecto. Se tendrá un código que corre sólo sobre CPU y otro que hace uso de GPU (CUDA). Luego de esto se debe hacer lo siguiente una vez se decida qué código se quiera compilar. En este caso se hará el ejemplo con el código en GPU.

    cd gpuSistemaieee118
    mkdir build
    cd build
    cmake ..
    make

Con estos sencillos pasos *cmake* buscará las dependencias y verá que se puedan cumplir, en caso negativo se le informará al usuario para que proceda con la instalación de las librerías necesarias.
