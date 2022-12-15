// gcc OWM.c -o func_OWM.o -lm -llas_c
//./datos.o ../datos/$(ls ../datos/*.las)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

// #include <liblas/capi/las_version.h> //porque esta se llama ya desde liblas.h
#include <liblas/capi/liblas.h>

// #include "utils/file.h"
// #include "utils/point.h"
// #include "utils/vector.h"
// #include "utils/util.h"
// #include "utils/plane.h"
// #include "utils/octree.h"
// #include "utils/main_options.h"

#include "include/environment.h"

// #define NUM_PROCS 4

// #pragma GCC push_options
// #pragma GCC optimize ("unroll-loops")

// __attribute__((optimize("unroll-loops")))

// unsigned int findMin(Lpoint** neighbors, double zzmin, unsigned int cellPoints) {
//   unsigned int idmin=0;
//   for(int i=1; i<cellPoints; i++)
//     if(neighbors[i]->z < zzmin){
//       zzmin = round(neighbors[i]->z * 100)/100;
//       idmin=i;
//     }
//   return idmin;
// }

// #pragma GCC pop_options

int main( int argc, char* argv[]){

    Vector3D max,min;

    LASReaderH reader=NULL;
    LASHeaderH header = NULL;

    unsigned int Npoints=0;

    // char inputTXT[50] = {"../datos/INAER_2011_Alcoy.xyz"};
    char inputLAS[50] = {"../datos/INAER_2011_Alcoy.las"};

    // Compruebo los argumentos de entrada

    for(int i=1 ; i < argc ; i++){
        strcpy(inputLAS,argv[i]);
        printf("Input.las: %s\n", inputLAS);
        reader = LASReader_Create(inputLAS);
        header = LASReader_GetHeader(reader);
        Npoints = LASHeader_GetPointRecordsCount(header);
        printf("Número de puntos      %d\n",Npoints);
        min.x = LASHeader_GetMinX(header), max.x = LASHeader_GetMaxX(header);
        min.y = LASHeader_GetMinY(header), max.y = LASHeader_GetMaxY(header);
        min.z = LASHeader_GetMinZ(header), max.z = LASHeader_GetMaxZ(header);
        printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );
    }
    // if(argc>1) {
    //   // strcpy(inputTXT,argv[1]);
    //   strcpy(inputLAS,argv[1]);
    //   // strcat(inputTXT,".xyz");
    //   // strcat(inputLAS,".las");
    // }

    // // omp_set_num_threads(num_procs);
    //
    // // printf("Input.txt: %s\n", inputTXT);
    // printf("Input.las: %s\n", inputLAS);
    //
    //
    // reader = LASReader_Create(inputLAS);
    // // reader = LASReader_Create("../datos/INAER_2011_Alcoy.las");
    // // reader = LASReader_Create("../datos/INAER_2011_Alcoy_Core.las");
    // // reader = LASReader_Create("../datos/INSITU_2018_Jenaro.las");
    // header = LASReader_GetHeader(reader);
    //
    //
    // // Recorro todo el file para contar las líneas
    // unsigned int Npoints=0;
    // Npoints = LASHeader_GetPointRecordsCount(header);
    // // while(!feof(fileLAS))
    // //   if(fgetc(fileLAS) == '\n')
    // //     Npoints++;
    // printf("Número de puntos      %d\n",Npoints);
    //
    //
    // min.x = LASHeader_GetMinX(header), max.x = LASHeader_GetMaxX(header);
    // min.y = LASHeader_GetMinY(header), max.y = LASHeader_GetMaxY(header);
    // min.z = LASHeader_GetMinZ(header), max.z = LASHeader_GetMaxZ(header);

    //Ya no necesito más el LAS
    LASReader_Destroy(reader);
    LASHeader_Destroy(header);




  return 0;
}
