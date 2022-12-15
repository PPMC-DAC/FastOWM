// gcc OWM.c -o func_OWM.o -lm -llas_c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

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



#define EPS 1e-10

int main( int argc, char* argv[]){


  LASReaderH reader=NULL;
  LASHeaderH header = NULL;

  // Tamaño de la ventana deslizante
  unsigned short Wsize = 12;
  // Tamaño de la rejilla
  unsigned char Bsize = 20;
  // Solape de la ventana deslizante
  float Overlap = 0;

  char inputTXT[50] = {"../datos/INSITU_2018_Jenaro.xyz"};
  char inputLAS[50] = {"../datos/INSITU_2018_Jenaro.las"};

  // Compruebo los argumentos de entrada
  // if(argc>1) Wsize = atoi(argv[0]);
  // if(argc>2) Bsize = atoi(argv[1]);
  // if(argc>3) Overlap = atoi(argv[2]);

  if(argc>1) {
    strcpy(inputTXT,argv[1]);
    strcpy(inputLAS,argv[1]);
    strcat(inputTXT,".xyz");
    strcat(inputLAS,".las");
  }

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
  Octree octreeInit = NULL;
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

  // for(int i=0; i<Npoints ; i++){
  //     if(getline(&line, &length, fileLAS) < 0){
  //       printf("File ended before all points could be read\n");
  //       return -1;
  //     }
  // }
  printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );

  radius = getRadius(min, max, &maxRadius);
  center = getCenter(min, radius);
  printf("Octree: Radio maximo con la mayor diferencia entre las dimensiones\n");
  printf("MaxRadius:  %.3f\n", maxRadius);
  printf("Center:     %.2f , %.2f\n", center.x,center.y);
  printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
  printf("CREANDO OCTREE...\n");
  octreeInit = createOctree(center, maxRadius);
  for(int i = 0; i < Npoints; i++)
	   // insertPointMinRadius(&pointer[i], octreeInit, main_options.octree_min_radius);
     insertPoint(&pointer[i], octreeInit);


  double Width = round((max.x-min.x)*100)/100;
  double High = round((max.y-min.y)*100)/100;
  double Density = Npoints/(Width*High);
  printf("Nube de puntos\n");
  printf("Ancho:  %.2lf\n",Width);
  printf("Alto:  %.2lf\n",High);
  printf("Densidad:  %.3lf\n",Density);


  unsigned short minNumPoints = 0.5*Density*Wsize*Wsize;
  printf("Minimo numero de puntos por celda:   %u\n", minNumPoints);

  unsigned short Displace = Wsize*(1-Overlap);
  printf("Displacement   %u\n", Displace);
  unsigned short NumDisplace = Wsize/Displace-1;
  printf("N displacements   %u\n", NumDisplace);

  unsigned short Crow=(int)round(Width/Wsize), Ccol=(int)round(High/Wsize);
  printf("Celdas por columa: %d\n",Ccol);
  printf("Celdas por fila:   %d\n",Crow);

  unsigned short Ncells=Crow*Ccol;
  printf("N cells:    %u\n", Ncells);

  clock_t t_init = clock();
  printf("/////////////////////////// LOOP ///////////////////////////\n\n\n");

  // int j,i;
  // Lcell myCell = {0.0,0.0,0.0,0.0,0};
  // for(j=-NumDisplace; j<=NumDisplace; j++ ){
  //   for(i=-NumDisplace; i<=NumDisplace; i++){
  //     printf("Point [%d, %d]\n",j,i );
  //     // myCell.xmin = pointer[Wsize];
  //     // myCell.xmax = pointer[0];
  //   }
  // }


  Lpoint** neighbors = NULL;
  // Octree cell = NULL;
  Lpoint cellCenter;
  double zzmin;
  unsigned int cellNum=0, cellPoints=0;
  unsigned int idmin = 0;
  cellCenter.y = min.y + Wsize/2;
  while(cellCenter.y < max.y){
      cellCenter.x = min.x + Wsize/2;
      while(cellCenter.x < max.x){
          //Este punto no lo compara con nada
          // *cellCenter = createPoint(0, cX, cY, 0.1, 0, zeroArray);
          // cellCenter->x = cX;
          // cellCenter->y = cY;
          // Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs)
          // printf("Busco los vecinos\n");
          neighbors = searchNeighbors2D(&cellCenter, octreeInit, Wsize/2, &cellPoints);
          // neighbors = searchOctNeighborsFelipe(&cellCenter, octreeInit, Wsize/2, &cellPoints);
          // printf("Numero de elementos de la celda: %d\n", cellPoints );
          //Tengo que hacerlo porque el algoritmo no lo hace
          if(cellPoints >= minNumPoints ){
              printf("Numero de elementos de la celda: %d\n", cellPoints );
              zzmin = neighbors[0]->z;
              // // printf("El mínimo de inicio es: %.2f\n", zzmin);
              // if(cellNum == 54)
              //   for(int i=0; i<cellPoints; i++){
              //       printf("Vecino: %.2f\n", neighbors[i]->z);
              //   }
              // center.x=cellCenter.x;
              // center.y=cellCenter.y;
              // cell = createOctree(center,Wsize/2);
              // for(int i = 0; i < cellPoints; i++)
            	//    insertPointMinRadius(neighbors[i], cell, 0.0);
              // searchMinFelipe(cell,&zzmin,&idmin);
              // // printf("El mínimo de la celda es: %.2f\n", pointer[idmin].z);
              idmin = 0;
              // while(abs(neighbors[idmin]->z - zzmin) > 0.0001 && idmin < cellPoints) idmin++;
              for(int i=1; i<cellPoints; i++)
                if(round(neighbors[i]->z * 100)/100 < zzmin){
                  // if(cellNum == 54) printf("Min potencial %.2f\n", neighbors[i]->z);
                  zzmin = neighbors[i]->z;
                  idmin=i;
                }
              // printf("El mínimo de la celda es: %.2f\n", neighbors[idmin]->z);
              cellNum++;
          }
          cellCenter.x += Wsize;
          free(neighbors);
          neighbors = NULL;
          // free(cell);
          // cell=NULL;
      }
      cellCenter.y += Wsize;
  }
  printf("\nCeldas no descartadas    %d\n", cellNum);
  // free(neighbors);
  // neighbors = NULL;
  // free(cell);
  // cell = NULL;

  printf("\n\n/////////////////////////// END ///////////////////////////\n");

  // Compruebo que mínimos se repiten
  if(NumDisplace > 0){

  }

  // Aplico la malla para ver si hay zonas sin puntos
  if(Bsize > 0){
    printf("Aplico la malla\n");
  }

  // Libero memoria y cierro el fichero
  free(pointer);
  free(octreeInit);
  if(fclose(fileLAS)){
    printf("Cannot close the file\n");
    return -1;
  }
  printf("Elapsed:     %.6f s\n",(double)(clock()-t_init) / CLOCKS_PER_SEC );
  return 0;
}
