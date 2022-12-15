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

    typedef struct{
        unsigned int id;
        unsigned int count;
    } accMin;


    double t_stage, t_func;

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

    unsigned short num_procs = 4;

    unsigned short block_size = 1;

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
    if(argc>6) block_size = atoi(argv[6]);

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

    // #pragma omp parallel for firstprivate(Npoints,pointer,octreeIn) schedule(static,1)
    for(int i = 0; i < Npoints; i++)
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

    // unsigned short Crow=(int)round(Width/Wsize), Ccol=(int)round(High/Wsize);
    unsigned short  Crow=(int)floor(Width/Wsize)+1, Ccol=(int)floor(High/Wsize)+1;
    if(Overlap!=0) Crow=(int)(round((Width+2*Wsize*Overlap)/Displace))-1, Ccol=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
    printf("Celdas por columa: %d\n",Ccol);
    printf("Celdas por fila:   %d\n",Crow);

    printf("\nN cells:    %u\n", Crow*Ccol);

    // Como voy a tener como mucho un mínimo por celda...
    int* minIDs = malloc(Crow*Ccol*sizeof(int));

    t_func=omp_get_wtime();
    t_stage=omp_get_wtime();
    printf("/////////////////////////// LOOP ///////////////////////////\n\n");

    Lpoint** neighbors = NULL;
    Lpoint cellCenter = {0,0.0,0.0,0.0};
    double zzmin;
    unsigned int index=0, countMin=0, cellPoints=0;
    unsigned int idmin = 0;

    /* El extremo izquierdo de la ventana es el que tiene que salirse de la nube!!
    Sino, me dejo puntos. Puedo utilizar un factor para controlar esto en vez de que
    siempre sea el mismo valor. Esto no tiene sentido si estoy haciendo el solapado,
    pero si tiene mucho si tengo Overlap = 0, para poder controlar la última posición de
    la ventana*/
    // Wfactor=1 : estoy considerando el extremo izquierdo de la ventana
    // Wfactor=0.5 : considero el centro
    // Wfactor=0 : el extremo derecho
    // double Wfactor = 0.8; // esto quiere decir que permito que el centro se salga un poco de la nube
    // if(Overlap!=0) Wfactor = 0;

    double initX = min.x + Wsize/2 - Wsize*Overlap;
    double initY = min.y + Wsize/2 - Wsize*Overlap;

    // int myid;
    // cellCenter.y = initY;
    // while(cellCenter.y + Wsize/2 - Wfactor*Wsize < max.y + Wsize*Overlap){
    int ii,jj;
    // int iteraciones=0;
    // #pragma omp parallel for shared(countMin) firstprivate(minIDs,octreeIn,Wsize,Overlap,Displace,initX,initY) \
    //                           private(cellCenter,neighbors,cellPoints,idmin) schedule(dynamic,block_size)
    #pragma omp parallel shared(countMin) firstprivate(Ccol,Crow,minIDs,octreeIn,Wsize,Displace,initX,initY) \
                              private(ii,jj,cellCenter,neighbors,cellPoints,idmin) num_threads(num_procs)
    {
        for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=num_procs ){
        // for( int ii = 0 ; cellCenter.y + Wsize/2 - Wfactor*Wsize < max.y + Wsize*Overlap ; cellCenter.y += Displace ){
            cellCenter.y = initY + jj*Displace;
            // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);
            // cellCenter.x = initX;
            // #pragma omp parallel shared(countMin,minIDs) firstprivate(Crow,cellCenter,octreeIn,Wsize,Overlap,Displace,initX) \
            //                           private(ii,neighbors,cellPoints,idmin) num_threads(1)
            // {
                // for( ii = omp_get_thread_num() ; ii < Crow ; ii+=1 ){
                #pragma omp parallel for shared(countMin) schedule(dynamic,block_size)
                for( ii=0 ; ii < Crow ; ii++ ){
                // while(cellCenter.x + Wsize/2 - Wfactor*Wsize < max.x + Wsize*Overlap){
                // for( cellCenter.x = initX ; cellCenter.x + Wsize/2 - Wfactor*Wsize < max.x + Wsize*Overlap ; cellCenter.x += Displace ){
                    cellCenter.x = initX + ii*Displace;
                    // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
                    // printf("Busco los vecinos\n");
                    neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
                    // printf("Numero de elementos de la celda: %d\n", cellPoints );
                    if(cellPoints >= minNumPoints ){
                        // printf("Numero de elementos de la celda: %d\n", cellPoints );
                        idmin = findMin(neighbors, round(neighbors[0]->z * 100)/100, cellPoints);
                        #pragma omp critical
                          countMin++;
                        // printf("El mínimo %d de la celda es: %.2f\n",countMin ,neighbors[idmin]->z);
                        minIDs[countMin-1] = neighbors[idmin]->id;
                    }
                    // cellCenter.x += Displace;
                    free(neighbors);
                    neighbors = NULL;
                }
            // }
            // cellCenter.y += Displace;
        }
    }

  printf("\nCeldas no descartadas:   %d\n", countMin);
  // printf("Minimos seleccionados:\n");
  // for(int i=0; i<countMin; i++) printf(" %.2f ", pointer[minIDs[i]].z); printf("\n");
  // for(int i=0; i<countMin; i++) printf(" %d ", minIDs[i]); printf("\n");

  printf("\n\n/////////////////////////// END ///////////////////////////\n");
  printf("Time elapsed at STAGE 1:    %.6f s\n\n", omp_get_wtime()-t_stage);

  // Lpoint** minimos = NULL;
  // Descarto mínimos
  accMin* vMin = malloc(countMin*sizeof(accMin));
  unsigned int searcher=countMin;
  double t_substage;
  // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
  if(Overlap != 0){
      t_stage=omp_get_wtime();
      printf("\nN minimos de entrada    %d\n", searcher);
      printf("/////////////////////////// MIN SELECT ///////////////////////////\n\n");

      for(int i=0; i<countMin; i++){
          //Compruebo dos cosas:
          //  1)Que no entcuentro la ID
          //  2)Que no he llegado al último de los almacenados
          // Si encuentro una ID, paro y acumulo. Si no tengo más guardados, paro y almaceno
          for( searcher=0 ; minIDs[i]!=vMin[searcher].id && searcher<index ; searcher++ );

          if(vMin[searcher].id == minIDs[i]) {
            vMin[searcher].count++;
            // printf("Acumulo el %d; Posicion %d\n", minIDs[i], searcher);
          }

          if(searcher==index) {
            // printf("NUEVO: %d; Posicion %d\n", minIDs[i], index);
            vMin[index].id = minIDs[i];
            vMin[index].count = 1;
            index++;
          }

          // searcher=0;
      }

      searcher=0;
      #pragma omp parallel for shared(searcher) schedule(dynamic,block_size)
      for(int i=0; i<index; i++){
        //Confío en que las posiciones que no se han modificado están a 0
        if(vMin[i].count > 1){
          #pragma omp critical
          {
            minIDs[searcher] = vMin[i].id;
            searcher++;
          }
        }
      }

      // // int i1;
      // // #pragma omp parallel private(i1,searcher) shared(vMin,minIDs) firstprivate(countMin)
      // // {
      //     // printf("thread: %d\n", omp_get_thread_num());
      //     // for(i1=omp_get_thread_num(); i1<countMin; i1+=num_procs){
      //     for(int i1=0; i1<countMin; i1++){
      //         // if(minIDs[i]==0) printf("\nHE ENCONTRADO UNA ID = 0\n");
      //         for( searcher=0 ; minIDs[i1]!=vMin[searcher].id && searcher<countMin ; searcher++ );
      //         //si me salgo del bucle porque he llegado al final -> NUEVO
      //         if(searcher==countMin) {
      //           // printf("NUEVO: %d; Posicion %d\n", minIDs[i], index);
      //           vMin[i1].id = minIDs[i1];
      //           vMin[i1].count = 1;
      //           index++;
      //         }
      //     }
      // // }
      //
      // #pragma omp parallel for shared(vMin) firstprivate(minIDs) private(searcher)
      // for(int i=0; i<countMin; i++){
      //     for( searcher=0 ; searcher<countMin ; searcher++ ){
      //         //si me salgo y las IDs son iguales -> incremento el contador
      //         if(vMin[i].id == minIDs[searcher]) {
      //           vMin[i].count++;
      //           // index++;
      //           // printf("Acumulo el %d; Posicion %d\n", minIDs[i], searcher);
      //         }
      //     }
      // }


      // // printf("Lo acumulado es: \n");
      // // for(int i=0; i<index; i++) printf(" %d ", vMin[i].count); printf("\n");
      //
      // searcher=0;
      //
      // for(int i=0; i<countMin; i++){
      //   //Confío en que las posiciones que no se han modificado están a 0
      //   if(vMin[i].count > 2){
      //     // #pragma omp critical
      //     // {
      //       minIDs[searcher] = vMin[i].id;
      //       searcher++;
      //     // }
      //   }
      // }

      // printf("Finalmente me quedo con: \n");
      // for(int i=0; i<countMin; i++) printf(" %.2f ", pointer[minIDs[i]].z); printf("\n");
      // for(int i=0; i<searcher ; i++) printf(" %d ", minIDs[i]); printf("\n");

      printf("\nNumero de minimos que me quedo: %d ; %d ; %d\n", searcher, index, countMin);

      free(vMin);
      vMin = NULL;
      printf("\n\n/////////////////////////// END ///////////////////////////\n");
      printf("Time elapsed at STAGE 2:     %.6f s\n\n",omp_get_wtime() - t_stage );
  }


  // Aplico la malla para ver si hay zonas sin puntos
  Octree grid = NULL;
  Lpoint** minimos = NULL;
  // unsigned int emptyCell=0, fillCell=0;
  unsigned int addMin=0;
  //Aquí es donde tengo que tomar la decisión de si voy con todo o desprecio el ultimo cacho de nube
  Crow=(int)floor(Width/Bsize)+1;
  Ccol=(int)floor(High/Bsize)+1;
  // Crow=(int)round(Width/Bsize);
  // Ccol=(int)round(High/Bsize);
  printf("Dimensiones de la malla %dx%d\n", Ccol,Crow);
  printf("N cells grid:    %u\n", Crow*Ccol);
  // int minGridIDs[Crow*Ccol];
  int* minGridIDs = malloc(Crow*Ccol*sizeof(int));
  if(Bsize > 0){
      t_stage=omp_get_wtime();
      printf("/////////////////////////// GRID ///////////////////////////\n\n\n");

      //Creo un nuevo octree con todos los mínimos; las mismas dimensiones que el grande
      radius = getRadius(min, max, &maxRadius);
      center = getCenter(min, radius);
      printf("Creo el octree\n");
      grid = createOctree(center, maxRadius);
      for(int i = 0; i < searcher; i++)
         insertPoint(&pointer[minIDs[i]], grid);
      printf("Termino de insertar los puntos\n");

      #pragma omp parallel shared(addMin) firstprivate(Ccol,Crow,minGridIDs,octreeIn,grid,Bsize) \
                                private(ii,jj,cellCenter,neighbors,minimos,cellPoints,idmin) num_threads(num_procs)
      {
          // while(cellCenter.y - 0*Bsize/2 < max.y){
          for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=num_procs ){
             cellCenter.y = min.y + Bsize/2 + jj*Bsize;
             // cellCenter.x = min.x + Bsize/2;
             // while(cellCenter.x - 0*Bsize/2 < max.x){
             for( ii = 0 ; ii < Crow ; ii++ ){
                 cellCenter.x = min.x + Bsize/2 + ii*Bsize;

                 minimos = searchNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
                 // printf("Numero de elementos de la celda: %d\n", cellPoints );
                 //Tengo que hacerlo porque el algoritmo no lo hace
                 if(cellPoints == 0){
                     // printf("    Tengo una celda vacia en la malla\n" );
                     neighbors = searchNeighbors2D(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                     if(cellPoints>0){

                         idmin = findMin(neighbors, round(neighbors[0]->z * 100)/100, cellPoints);

                         minGridIDs[addMin] = neighbors[idmin]->id;
                         addMin++;
                         // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
                     }
                     // emptyCell++;
                     free(neighbors);
                     neighbors = NULL;

                 }//else fillCell++;
                 // cellCenter.x += Bsize;
                 free(minimos);
                 minimos=NULL;
             }
             // cellCenter.y += Bsize;
          }
      }
      printf("Minimos añadidos:        %d\n", addMin);
      // for(int i=0; i<addMin ; i++)
      //   printf("%d: %.2f , %.2f , %.2f\n", i+1, pointer[minGridIDs[i]].x, pointer[minGridIDs[i]].y,pointer[minGridIDs[i]].z);
      printf("\n\n/////////////////////////// END ///////////////////////////\n");
      printf("Time elapsed at STAGE 3     %.6f s\n\n",omp_get_wtime() - t_stage );
      free(grid);
      grid = NULL;
  }
  // printf("Minimos añad.    %d\n", addMin);
  // printf("Celdas vacias    %d\n", emptyCell);
  // printf("Celdas llenas    %d\n", fillCell);
  // printf("Total            %d: OK? %d\n", emptyCell+fillCell, emptyCell+fillCell == Crow*Ccol);

  printf("Finalmente, el mapa va a tener %d puntos, %d puntos menos\n", searcher+addMin, Npoints - searcher+addMin );
  // Ya no necesito el octreeIn
  free(octreeIn);
  octreeIn=NULL;

  FILE* fileMin;
  char outputTXT[50] = {"../datos/INAER_2011_Alcoy_salida.xyz"};
  char outputLAS[50] = {"../datos/INAER_2011_Alcoy_salida.LAS"};
  // Compruebo los argumentos de entrada
  if(argc>1) {
    strcpy(outputTXT,argv[1]);
    strcpy(outputLAS,argv[1]);
    strcat(outputTXT,"_salida.xyz");
    strcat(outputLAS,"_salida.las");
  }

  printf("Creo el fichero %s ...\n", outputTXT);
  if((fileMin = fopen(outputTXT,"w")) == NULL){
    printf("Unable to create file!\n");
    return -1;
  }

  for(int i=0; i<searcher;i++)
    fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);

  for(int i=0;i<addMin;i++)
    fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minGridIDs[i]].x, pointer[minGridIDs[i]].y,pointer[minGridIDs[i]].z);

  //Ya no necesito mas el fichero
  if(fclose(fileMin)){
    printf("Cannot close the file\n");
    return -1;
  }

  printf("Total time elapsed:     %.6f s\n", omp_get_wtime() - t_func);

  // Libero memoria
  free(minIDs);
  free(minGridIDs);
  free(pointer);


  // char command[200];
  printf("Creo el fichero %s ...\n", outputLAS);
  printf("txt2las -i %s -o %s --parse xyz\n",outputTXT,outputLAS);
  printf("> source(\"/home/felipe/Escritorio/Beca_CiTIUS/Proyecto_LiDAR/algoritmo_OWM_LiDAR/pintaLAS.R.R\")\n");
  // sprintf(command,"txt2las -i %s -o %s --parse xyz",outputTXT,outputLAS);
  // // sprintf(command,"python3 script_surface.py %s %s",outputTXT,outputLAS);
  // if(system(command) < 0){
  //   perror("system");
  //   return 1;
  // }


  return 0;
}
