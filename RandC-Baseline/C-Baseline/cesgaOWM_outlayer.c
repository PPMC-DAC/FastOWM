
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

#include "include/environment.h"

int cmpLpoint (const void* a, const void* b){
   // Le pongo el 100 para que se quede con el signo
   // return 100*( round( ((zSearch*)a)->z * 100)/100 - round( ((zSearch*)b)->z * 100)/100 );
   Lpoint* punteroA = (Lpoint*) a;
   Lpoint* punteroB = (Lpoint*) b;

   if (round2d( (punteroB)->z ) > round2d( (punteroA)->z )) {
        return -1;
    } else if (round2d( (punteroB)->z ) < round2d( (punteroA)->z )) {
        return 1;
    } else {
        return 0;
    }

}

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int main( int argc, char* argv[]){

    // typedef struct{
    //     unsigned int id;
    //     unsigned short count;
    // } accMin;

    // Ficheros
    FILE* fileLAS;
    FILE* fileMin;

    // Octrees
    Octree octreeIn = NULL;
    Octree grid = NULL;

    //Listas
    Lpoint* pointer = NULL;
    int* minIDs = NULL;
    // accMin* vMin = NULL;
    int* minGridIDs = NULL;
    Lpoint** minimos = NULL;

    unsigned int Npoints=0, Ncells, Ngrid;
    Vector3D center, radius, min, max;
    float maxRadius = 0.0;
    double Width, High, Density, Displace;
    unsigned short  minNumPoints, Crow, Ccol, Crowg, Ccolg;

    Lpoint cellCenter = {0,0.0,0.0,0.0};
    // double zzmin;
    unsigned int index=0, countMin=0, cellPoints=0;
    // unsigned int idmin = 0;

    // double initX, initY;
    int ii,jj; //para el control de los bucles paralelizados

    // unsigned int countMin;
    //
    // unsigned int searcher;

    unsigned int addMin=0;

    double t_stage, t_func;

    // Tamaño de la ventana deslizante
    // unsigned short Wsize = 12;
    // Tamaño de la rejilla
    unsigned short Gsize = 20;
    // Factor para eliminar el minimo
    double out_factor = 0.8;
    // Numero de procesadores
    // unsigned short num_procs = 4;
    // Tamaño del bloque del scheduler de OMP
    // unsigned short block_size = 1;

    //Control del bucle de ejecución y parámetros
    // unsigned int bucle_entrada = 1;

    char inputTXT[50] = {"../datos/INAER_2011_Alcoy.xyz"};
    // char inputLAS[50] = {"../datos/INAER_2011_Alcoy.las"};
    char outputTXT[50] = {"../datos/INAER_2011_Alcoy_salida.xyz"};
    // char outputLAS[50] = {"../datos/INAER_2011_Alcoy_salida.las"};

    // Compruebo los argumentos de entrada
    if(argc>1) {
      strcpy(inputTXT,argv[1]);
      // strcpy(inputLAS,argv[1]);
      strcat(inputTXT,".xyz");
      // strcat(inputLAS,".las");

      strcpy(outputTXT,argv[1]);
      // strcpy(outputLAS,argv[1]);
      strcat(outputTXT,"_clean.xyz");
      // strcat(outputLAS,"_salida.las");
    }
    if(argc>2) Gsize = atoi(argv[2]);
    // if(argc>3) Bsize = atoi(argv[3]);
    if(argc>3) out_factor = atof(argv[3]);
    // if(argc>3) num_procs = atoi(argv[3]);
    // if(argc>6) bucle_entrada = atoi(argv[6]);

    printf("Tamaño de rejilla     %u\n", Gsize);

    // resultados = malloc(bucle_entrada*sizeof(double));
    // double resultados[bucle_entrada];

    // omp_set_num_threads(num_procs);

    printf("Input.txt: %s \n", inputTXT);
    // printf("Input.las: %s\n", inputLAS);

    // Abro el fichero
    if((fileLAS = fopen(inputTXT,"r")) == NULL){
      printf("Unable to open file!\n");
      return -1;
    }

    // // Numero de puntos para hacer la reserva de memoria justa
    while(!feof(fileLAS))
      if(fgetc(fileLAS) == '\n')
        Npoints++;
    // Vuelvo al principio del fichero
    rewind(fileLAS);
    printf("Número de puntos      %d\n",Npoints);

    // Reservo memoria para la nube de puntos
    pointer = malloc(Npoints*sizeof(Lpoint));

    printf("VOLCANDO PUNTOS...\n");
    //Obtengo los datos id X Y Z
    pointer[0].id = 0;
    fscanf(fileLAS, "%lf%lf%lf",&pointer[0].x,&pointer[0].y,&pointer[0].z);
    // Los inicializo así para que no haya problemas con los mínimos
    min.x=pointer[0].x;
    min.y=pointer[0].y;
    //Los siguientes 6 campos no los necesito
    while(fgetc(fileLAS)!='\n');
    // printf("[%u]  %.2lf, %.2lf, %.2lf\n",pointer[0].id,pointer[0].x,pointer[0].y,pointer[0].z);

    for(int i=1; i<Npoints ; i++){
      //Obtengo los datos id X Y Z
      pointer[i].id = i;
      if(fscanf(fileLAS, "%lf%lf%lf",&pointer[i].x,&pointer[i].y,&pointer[i].z) < 3){
        printf("Imposible to obtain values\n");
        return -1;
      }
      if(min.x>pointer[i].x) min.x=pointer[i].x;
      if(max.x<pointer[i].x) max.x=pointer[i].x;
      if(min.y>pointer[i].y) min.y=pointer[i].y;
      if(max.y<pointer[i].y) max.y=pointer[i].y;
      //Los siguientes 6 campos no los necesito
      while(fgetc(fileLAS)!='\n');
      // printf("[%u]  %.2lf, %.2lf, %.2lf\n",pointer[i].id,pointer[i].x,pointer[i].y,pointer[i].z);
    }

    //Ya no necesito mas el fichero
    if(fclose(fileLAS)){
      printf("Cannot close the file\n");
      return -1;
    }

    // qsort(pointer,Npoints,sizeof(Lpoint),&cmpLpoint);
    // for(int i=0; i<Npoints; i++) printf(" %.2f ", pointer[i].z); printf("\n\n");

    Width = round2d(max.x-min.x);
    High = round2d(max.y-min.y);
    printf("Ancho:  %.2lf\n",Width);
    printf("Alto:  %.2lf\n",High);

    //Aquí es donde tengo que tomar la decisión de si voy con todo o desprecio el ultimo cacho de nube
    Crow=(int)floor(Width/Gsize)+1;
    Ccol=(int)floor(High/Gsize)+1;
    // Crow=(int)round(Width/Bsize);
    // Ccol=(int)round(High/Bsize);
    printf("Dimensiones de la malla %dx%d\n", Ccol,Crow);
    printf("N cells grid:    %u\n", Crow*Ccol);

    minIDs = malloc(Npoints*sizeof(int));

    t_stage=omp_get_wtime();
    printf("/////////////////////////// GRID ///////////////////////////\n\n\n");

    //Creo un nuevo octree con todos los mínimos; las mismas dimensiones que el grande
    radius = getRadius(min, max, &maxRadius);
    center = getCenter(min, radius);
    // printf("Creo el octree\n");
    grid = createOctree(center, maxRadius);
    for(int i = 0; i < Npoints; i++)
       insertPoint(&pointer[i], grid);
    // printf("Termino de insertar los puntos\n");

    // addMin = stage3(Bsize, Crow, Ccol, minGridIDs, octreeIn, grid, min);

    double umbral = 0.0;
    // #pragma omp parallel shared(addMin) firstprivate(Ccol,Crow,grid,Gsize) \
    //                           private(ii,jj,cellCenter,neighbors,minimos,cellPoints,idmin)
    // {

        // for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=num_procs ){
        for( jj = 0 ; jj < Ccol ; jj++ ){
           cellCenter.y = min.y + Gsize/2 + jj*Gsize;

           for( ii = 0 ; ii < Crow ; ii++ ){
               umbral = 0.0;
               cellCenter.x = min.x + Gsize/2 + ii*Gsize;

               minimos = searchNeighbors2D(&cellCenter, grid, Gsize/2, &cellPoints);
               // printf("Numero de elementos de la celda: %d\n", cellPoints );
               // Hago la media de todas las alturas de la zona
               if(cellPoints > 0){
                   for( int i=0; i<cellPoints ; i++){
                     umbral += minimos[i]->z;
                   }
                   umbral = round2d(out_factor * umbral/cellPoints);
                   // printf("Umbral: %.2f\n", umbral);
                   // Voy comprobando cual es inferior a este
                   for( int i=0 ; i<cellPoints ; i++){
                     // Si es inferior lo acumulo en el 0 (porque lo he decidido así)
                     if(round2d(minimos[i]->z) > umbral){
                       minIDs[addMin] = minimos[i]->id;
                       addMin++;
                     }
                   }
                }

               free(minimos);
               minimos=NULL;
           }
        }
    // }
    printf("Minimos con los que me quedo:        %d\n", addMin);
    // for(int i=0; i<addMin ; i++)
    //   printf("%d: %.2f , %.2f , %.2f\n", i+1, pointer[minGridIDs[i]].x, pointer[minGridIDs[i]].y,pointer[minGridIDs[i]].z);
    printf("\n\n/////////////////////////// END ///////////////////////////\n");
    printf("Time elapsed at STAGE 3:     %.6f s\n\n",omp_get_wtime() - t_stage );
    // Ya no necesito este octree
    free(grid);
    grid = NULL;

    // Fichero de salida
    printf("Creo el fichero %s ...\n", outputTXT);
    if((fileMin = fopen(outputTXT,"w")) == NULL){
      printf("Unable to create file!\n");
      return -1;
    }

    for(int i=0 ; i<addMin ; i++)
      fprintf(fileMin, "%.2lf %.2lf %.2lf\n", pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);

    //Ya no necesito mas el fichero
    if(fclose(fileMin)){
      printf("Cannot close the file\n");
      return -1;
    }


    return 0;
}
