// gcc OWM.c -o func_OWM.o -lm -llas_c

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

#define NUM_PROCS 4

// #pragma GCC push_options
// #pragma GCC optimize ("unroll-loops")

// __attribute__((optimize("unroll-loops")))

unsigned int findMin(Lpoint** neighbors, double zzmin, unsigned int cellPoints) {
  unsigned int idmin=0;
  for(int i=1; i<cellPoints; i++)
    if(neighbors[i]->z < zzmin){
      zzmin = round(neighbors[i]->z * 100)/100;
      idmin=i;
    }
  return idmin;
}

// #pragma GCC pop_options

int main( int argc, char* argv[]){

    typedef struct{
        unsigned int id;
        unsigned int count;
    } accMin;


    double t1, t2, t3;

    LASReaderH reader=NULL;
    LASHeaderH header = NULL;

    // Tamaño de la ventana deslizante
    unsigned short Wsize = 12;
    // Tamaño de la rejilla
    unsigned short Bsize = 20;
    // Solape de la ventana deslizante
    double Overlap = 0.8;

    char inputTXT[50] = {"../datos/INAER_2011_Alcoy.xyz"};
    char inputLAS[50] = {"../datos/INAER_2011_Alcoy.las"};

    unsigned short num_procs = 5;

    // Compruebo los argumentos de entrada
    if(argc>1) {
      strcpy(inputTXT,argv[1]);
      strcpy(inputLAS,argv[1]);
      strcat(inputTXT,".xyz");
      strcat(inputLAS,".las");
    }
    if(argc>2) Wsize = atoi(argv[2]);
    if(argc>3) Bsize = atoi(argv[3]);
    if(argc>4) Overlap = atof(argv[4]);
    if(argc>5) num_procs = atoi(argv[5]);

    omp_set_num_threads(num_procs);

    printf("Input.txt: %s\n", inputTXT);
    printf("Input.las: %s\n", inputLAS);


    reader = LASReader_Create(inputLAS);
    // reader = LASReader_Create("../datos/INAER_2011_Alcoy.las");
    // reader = LASReader_Create("../datos/INAER_2011_Alcoy_Core.las");
    // reader = LASReader_Create("../datos/INSITU_2018_Jenaro.las");
    header = LASReader_GetHeader(reader);

    // Abro el fichero
    FILE* fileLAS;
    if((fileLAS = fopen(inputTXT,"r")) == NULL){
    // if((fileLAS = fopen("../datos/INAER_2011_Alcoy.xyz","r")) == NULL){
    // if((fileLAS = fopen("../datos/INAER_2011_Alcoy_Core.xyz","r")) == NULL){
    // if((fileLAS = fopen("../datos/INSITU_2018_Jenaro.xyz","r")) == NULL){
      printf("Unable to open file!\n");
      return -1;
    }



    printf("PARAMETERS:\nTamaño de ventana     %u\n", Wsize);
    printf("Tamaño de rejilla     %u\nSolape                %.2f\n", Bsize,Overlap);

    // Recorro todo el file para contar las líneas
    unsigned int Npoints=0;
    Npoints = LASHeader_GetPointRecordsCount(header);
    // while(!feof(fileLAS))
    //   if(fgetc(fileLAS) == '\n')
    //     Npoints++;
    printf("Número de puntos      %d\n",Npoints);
    //
    // // Vuelvo al principio del fichero
    // rewind(fileLAS);

    Vector3D center, min, max, radius;
    float maxRadius = 0.0;
    double zeroArray[3] = {0.0,0.0,0.0};
    Octree octreeIn = NULL;
    // Reservo memoria para todos los puntos y obtengo los datos junto con los max y min
    // double xmax=LASHeader_GetMaxX(header),
    //        xmin=LASHeader_GetMinX(header),
    //        ymax=LASHeader_GetMaxY(header),
    //        ymin=LASHeader_GetMinY(header);
    min.x = LASHeader_GetMinX(header), max.x = LASHeader_GetMaxX(header);
    min.y = LASHeader_GetMinY(header), max.y = LASHeader_GetMaxY(header);
    min.z = LASHeader_GetMinZ(header), max.z = LASHeader_GetMaxZ(header);

    //Ya no necesito más el LAS
    LASReader_Destroy(reader);
    LASHeader_Destroy(header);
    // Lpoint xmin, xmax, ymax, ymin;
    Lpoint* pointer = malloc(Npoints*sizeof(Lpoint));

    //Obtengo los datos id X Y Z
    // pointer[0].id = 0;
    // fscanf(fileLAS, "%lf%lf%lf",&pointer[0].x,&pointer[0].y,&pointer[0].z);
    // Los inicializo así para que no haya problemas con los mínimos
    // xmin=pointer[0].x;
    // xmax=pointer[0].x;
    // ymin=pointer[0].y;
    // ymax=pointer[0].y;
    //Los siguientes 6 campos no los necesito
    // while(fgetc(fileLAS)!='\n');

    printf("VOLCANDO PUNTOS...\n");
    for(int i=0; i<Npoints ; i++){
      //Obtengo los datos id X Y Z
      // pointer[i] = createPoint(i,0.0, 0.0, 0.0);
      pointer[i].id = i;
      if(fscanf(fileLAS, "%lf%lf%lf",&pointer[i].x,&pointer[i].y,&pointer[i].z) < 3){
        printf("Imposible to obtain values\n");
        return -1;
      }
      // if(xmin>pointer[i].x) xmin=pointer[i].x;
      // if(xmax<pointer[i].x) xmax=pointer[i].x;
      // if(ymin>pointer[i].y) ymin=pointer[i].y;
      // if(ymax<pointer[i].y) ymax=pointer[i].y;
      //Los siguientes 6 campos no los necesito
      while(fgetc(fileLAS)!='\n');
     	// printf("[%u]  %.2lf, %.2lf, %.2lf\n",pointer[i].id,pointer[i].x,pointer[i].y,pointer[i].z);
    }

    //Ya no necesito mas el fichero
    if(fclose(fileLAS)){
      printf("Cannot close the file\n");
      return -1;
    }

    printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );

    radius = getRadius(min, max, &maxRadius);
    center = getCenter(min, radius);
    printf("Octree: Radio maximo con la mayor diferencia entre las dimensiones\n");
    printf("MaxRadius:  %.3f\n", maxRadius);
    printf("Center:     %.2f , %.2f\n", center.x,center.y);
    printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
    printf("CREANDO OCTREE...\n");
    octreeIn = createOctree(center, maxRadius);

    for(int i = 0; i < Npoints; i++)
  	   // insertPointMinRadius(&pointer[i], octreeIn, main_options.octree_min_radius);
       insertPoint(&pointer[i], octreeIn);



    double Width = round((max.x-min.x)*100)/100;
    double High = round((max.y-min.y)*100)/100;
    double Density = Npoints/(Width*High);
    printf("Nube de puntos\n");
    printf("Ancho:  %.2lf\n",Width);
    printf("Alto:  %.2lf\n",High);
    printf("Densidad:  %.3lf\n",Density);


    unsigned short minNumPoints = 0.5*Density*Wsize*Wsize;
    printf("Minimo numero de puntos por celda:   %u\n", minNumPoints);

    float Displace = round(Wsize*(1-Overlap)*100)/100;
    printf("Displacement   %.2f\n", Displace);
    // unsigned short NumDisplace = Wsize/Displace-1;
    // printf("N displacements   %u\n", NumDisplace);

    unsigned short Crow=(int)round(Width/Wsize), Ccol=(int)round(High/Wsize);
    if(Overlap!=0) Crow=(int)(round((Width+2*Wsize*Overlap)/Displace))-1, Ccol=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
    printf("Celdas por columa: %d\n",Ccol);
    printf("Celdas por fila:   %d\n",Crow);

    printf("\nN cells:    %u\n", Crow*Ccol);

    // Como voy a tener como mucho un mínimo por celda...
    int* minIDs = malloc(Crow*Ccol*sizeof(int));

    clock_t t_phase, t_init = clock();
    t1=omp_get_wtime();
    printf("/////////////////////////// LOOP ///////////////////////////\n\n");

    Lpoint** neighbors = NULL;
    Lpoint cellCenter = {0,0.0,0.0,0.0};
    double zzmin;
    unsigned int countMin=0, cellPoints=0;
    unsigned int idmin = 0;

    /* El extremo izquierdo de la ventana es el que tiene que salirse de la nube!!
    Sino, me dejo puntos. Puedo utilizar un factor para controlar esto en vez de que
    siempre sea el mismo valor. Esto no tiene sentido si estoy haciendo el solapado,
    pero si tiene mucho si tengo Overlap = 0, para poder controlar la última posición de
    la ventana*/
    // Wfactor=1 : estoy considerando el extremo izquierdo de la ventana
    // Wfactor=0.5 : considero el centro
    // Wfactor=0 : el extremo derecho
    double Wfactor = 0.8; // esto quiere decir que permito que el centro se salga un poco de la nube
    if(Overlap!=0) Wfactor = 0;

    double initX = min.x + Wsize/2 - Wsize*Overlap;
    double initY = min.y + Wsize/2 - Wsize*Overlap;

    int myid;
    // cellCenter.y = initY;
    // while(cellCenter.y + Wsize/2 - Wfactor*Wsize < max.y + Wsize*Overlap){
    #pragma omp parallel for shared(countMin) firstprivate(minIDs,octreeIn,Wfactor,Wsize,Overlap,Displace,initX,initY) \
                              private(cellCenter,neighbors,cellPoints) schedule(dynamic,1)
    for( int ii = 0 ; ii < Ccol ; ii++ ){
    // for( int ii = 0 ; cellCenter.y + Wsize/2 - Wfactor*Wsize < max.y + Wsize*Overlap ; cellCenter.y += Displace ){
        cellCenter.y = initY + ii*Displace;
        // cellCenter.x = initX;
        // while(cellCenter.x + Wsize/2 - Wfactor*Wsize < max.x + Wsize*Overlap){
        for( cellCenter.x = initX ; cellCenter.x + Wsize/2 - Wfactor*Wsize < max.x + Wsize*Overlap ; cellCenter.x += Displace ){
            // printf("Busco los vecinos\n");
            // #pragma omp critical
              neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
            // printf("Numero de elementos de la celda: %d\n", cellPoints );
            if(cellPoints >= minNumPoints ){
                // printf("Numero de elementos de la celda: %d\n", cellPoints );
                idmin = findMin(neighbors, round(neighbors[0]->z * 100)/100, cellPoints);
                countMin++;
                // printf("El mínimo %d de la celda es: %.2f\n",countMin ,neighbors[idmin]->z);
                minIDs[countMin-1] = neighbors[idmin]->id;
            }
            // cellCenter.x += Displace;
            free(neighbors);
            neighbors = NULL;
        }
        // cellCenter.y += Displace;
    }

  printf("\nCeldas no descartadas:   %d\n", countMin);
  // printf("Minimos seleccionados:\n");
  // for(int i=0; i<countMin; i++) printf(" %.2f ", pointer[minIDs[i]].z); printf("\n");
  // for(int i=0; i<countMin; i++) printf(" %d ", minIDs[i]); printf("\n");

  printf("\n\n/////////////////////////// END ///////////////////////////\n");
  t2=omp_get_wtime();
  printf("Elpased:    %.6f s\n\n", t2-t1);

  // Lpoint** minimos = NULL;
  // Descarto mínimos
  accMin* vMin = malloc(countMin*sizeof(accMin));
  unsigned int index=0, searcher=countMin;
  // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
  if(Overlap != 0){
      t_phase=clock();
      printf("\nN minimos de entrada    %d\n", countMin);
      printf("/////////////////////////// MIN SELECT ///////////////////////////\n\n");
      searcher=0;
      for(int i=0; i<countMin; i++){
          //Compruebo dos cosas:
          //  1)Que no entcuentro la ID
          //  2)Que no he llegado al último de los almacenados
          // Si encuentro una ID, paro y acumulo. Si no tengo más guardados, paro y almaceno
          while(minIDs[i]!=vMin[searcher].id && searcher<index) searcher++;
          if(minIDs[i]==vMin[searcher].id) {
            vMin[searcher].count++;
            // printf("Acumulo el %d; Posicion %d\n", minIDs[i], searcher);
          }
          if(searcher==index) {
            // printf("NUEVO: %d; Posicion %d\n", minIDs[i], index);
            vMin[index].id = minIDs[i];
            vMin[index].count = 1;
            index++;
          }
          searcher=0;
      }
      // printf("Lo acumulado es: \n");
      // for(int i=0; i<index; i++) printf(" %d ", vMin[i].count); printf("\n");
      // Descarto todos aquellos que se han repetido solo una vez
      // Borro todo los que tenía
      free(minIDs);
      searcher=0;
      // Y añado solo los que necesito; Como maximo todos los que se repiten al menos una vez
      minIDs = malloc(index*sizeof(int));
      for(int i=0; i<index; i++){
        if(vMin[i].count > 1){
          minIDs[searcher] = vMin[i].id;
          searcher++;
        }
      }
      // printf("Finalmente me quedo con: \n");
      // for(int i=0; i<countMin; i++) printf(" %.2f ", pointer[minIDs[i]].z); printf("\n");
      // for(int i=0; i<searcher ; i++) printf(" %d ", minIDs[i]); printf("\n");

      printf("\nNumero de minimos que me quedo: %d\n", searcher);

      free(vMin);
      vMin = NULL;
      printf("\n\n/////////////////////////// END ///////////////////////////\n");
      printf("Elapsed:     %.6f s\n\n",(double)(clock()-t_phase) / CLOCKS_PER_SEC );
  }


  // Aplico la malla para ver si hay zonas sin puntos
  Octree grid = NULL;
  Lpoint** minimos = NULL;
  // unsigned int emptyCell=0, fillCell=0;
  unsigned int addMin=0;
  if(Bsize > 0){
      t_phase=clock();
      Crow=(int)floor(Width/Bsize)+1;
      Ccol=(int)floor(High/Bsize)+1;
      printf("N cells grid:    %u\n", Crow*Ccol);
      // int minGridIDs[Crow*Ccol];
      int* minGridIDs = malloc(Crow*Ccol*sizeof(int));
      printf("/////////////////////////// GRID ///////////////////////////\n\n\n");

      //Creo un nuevo octree con todos los mínimos; las mismas dimensiones que el grande
      radius = getRadius(min, max, &maxRadius);
      center = getCenter(min, radius);
      // printf("Creo el octree\n");
      grid = createOctree(center, maxRadius);
      for(int i = 0; i < searcher; i++)
         insertPoint(&pointer[minIDs[i]], grid);
      // printf("Termino de insertar los puntos\n");

      cellCenter.y = min.y + Bsize/2;
      // Voy a suponer que cuando el centro pase cierto valor el máximo de la ventana, ya la puedo descartar
      // Si es 0 quiere decir que considero que si el centro de la celda sobrepasa el maximo de la nube,
      // no tengo en cuenta esa ventana
      while(cellCenter.y - 0*Bsize/2 < max.y){
         cellCenter.x = min.x + Bsize/2;
         while(cellCenter.x - 0*Bsize/2 < max.x){

             minimos = searchNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
             // printf("Numero de elementos de la celda: %d\n", cellPoints );
             //Tengo que hacerlo porque el algoritmo no lo hace
             if(cellPoints == 0){
                 // printf("    Tengo una celda vacia en la malla\n" );
                 neighbors = searchNeighbors2D(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                 if(cellPoints>0){
                   // printf("    Tengo una celda vacia en la malla\n" );

                     // zzmin = round(neighbors[0]->z * 100)/100;
                     // idmin=0;
                     /* Imprescindible cuando hay muchos puntos.
                     Si da la casualidad de que es el primero y no he puesto idmin=0
                     si el numero de elementos de la celda es inferior al indice del minimo
                     de la celda anterior, esto PETA!! */
                     // for(int i=1; i<cellPoints; i++){
                     //      // printf("El minimo %.2f\n", zzmin);
                     //     if(neighbors[i]->z < zzmin){
                     //       zzmin = round(neighbors[i]->z * 100)/100;
                     //       idmin=i;
                     //     }
                     //   }

                     idmin = findMin(neighbors, round(neighbors[0]->z * 100)/100, cellPoints);

                     minGridIDs[addMin] = neighbors[idmin]->id;
                     addMin++;
                     // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
                 }
                 // emptyCell++;
                 free(neighbors);
                 neighbors = NULL;

             }//else fillCell++;
             cellCenter.x += Bsize;
             free(minimos);
             minimos=NULL;
         }
         cellCenter.y += Bsize;
      }
      printf("Minimos añadidos:        %d\n", addMin);
      // for(int i=0; i<addMin ; i++)
      //   printf("%d: %.2f , %.2f , %.2f\n", i+1, pointer[minGridIDs[i]].x, pointer[minGridIDs[i]].y,pointer[minGridIDs[i]].z);
      printf("\n\n/////////////////////////// END ///////////////////////////\n");
      printf("Elapsed:     %.6f s\n\n",(double)(clock()-t_phase) / CLOCKS_PER_SEC );
  }
  // printf("Minimos añad.    %d\n", addMin);
  // printf("Celdas vacias    %d\n", emptyCell);
  // printf("Celdas llenas    %d\n", fillCell);
  // printf("Total            %d: OK? %d\n", emptyCell+fillCell, emptyCell+fillCell == Crow*Ccol);

  printf("Finalmente, el mapa va a tener %d puntos, %d puntos menos\n", searcher+addMin, Npoints - searcher+addMin );


  // Libero memoria
  free(pointer);
  free(octreeIn);
  free(grid);

  printf("Total elapsed:     %.6f s\n",(double)(clock()-t_init) / CLOCKS_PER_SEC );
  return 0;
}
