// gcc OWM.c -o func_OWM.o -lm -llas_c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <liblas/capi/las_version.h> //porque esta se llama ya desde liblas.h
#include <liblas/capi/liblas.h>

#define EPS 1e-4

int main( int argc, char* argv[]){

  typedef struct{
    unsigned int id;
    double x;
    double y;
    double z;
  } Lpoint;

  typedef struct{
    unsigned int counter;
    Lpoint *pointer;
  } Lcloud;

  typedef struct{
    Lpoint xmin;
    Lpoint xmax;
    Lpoint ymin;
    Lpoint ymax;
    unsigned short Npoints;
  } Lcell;

  LASReaderH reader=NULL;
  LASHeaderH header = NULL;
  LASPointH p = NULL;
  reader = LASReader_Create("../datos/sample24.las");
  header = LASReader_GetHeader(reader);
  p = LASReader_GetNextPoint(reader);

  // Tamaño de la ventana deslizante
  unsigned short Wsize = 12;
  // Tamaño de la rejilla
  unsigned char Bsize = 20;
  // Solape de la ventana deslizante
  float Overlap = 0;

  clock_t t_init = clock();


  // Compruebo los argumentos de entrada
  if(argc>1) Wsize = atoi(argv[1]);
  if(argc>2) Bsize = atoi(argv[2]);
  if(argc>3) Overlap = atoi(argv[3]);

  printf("PARAMETERS:\nTamaño de ventana     %u\n", Wsize);
  printf("Tamaño de rejilla     %u\nSolape                %.2f\n", Bsize,Overlap);

  // Abro el fichero
  // FILE* fileLAS;
  // if((fileLAS = fopen("../datos/sample24.xyz","r")) == NULL){
  // // if((fileLAS = fopen("../datos/INAER_2011_Alcoy.xyz","r")) == NULL){
  //   printf("Unable to open file!\n");
	// 	return -1;
  // }

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

  // Reservo memoria para todos los puntos y obtengo los datos junto con los max y min
  double xmax=LASHeader_GetMaxX(header),
         xmin=LASHeader_GetMinX(header),
         ymax=LASHeader_GetMaxY(header),
         ymin=LASHeader_GetMinY(header);
  // Lpoint xmin, xmax, ymax, ymin;
  // Lcloud* cloud = malloc(Npoints* sizeof(Lcloud));

  //Obtengo los datos id X Y Z
  // pointer[0].id = 0;
  // fscanf(fileLAS, "%lf%lf%lf",&pointer[0].x,&pointer[0].y,&pointer[0].z);
  // // Los inicializo así para que no haya problemas con los mínimos
  // xmin=pointer[0].x;
  // xmax=pointer[0].x;
  // ymin=pointer[0].y;
  // ymax=pointer[0].y;
  // //Los siguientes 6 campos no los necesito
  // while(fgetc(fileLAS)!='\n');
  //
  // for(int i=1; i<Npoints-1 ; i++){
  //   //Obtengo los datos id X Y Z
  //   pointer[i].id = i;
  //   fscanf(fileLAS, "%lf%lf%lf",&pointer[i].x,&pointer[i].y,&pointer[i].z);
  //   if(xmin>pointer[i].x) xmin=pointer[i].x;
  //   if(xmax<pointer[i].x) xmax=pointer[i].x;
  //   if(ymin>pointer[i].y) ymin=pointer[i].y;
  //   if(ymax<pointer[i].y) ymax=pointer[i].y;
  //   //Los siguientes 6 campos no los necesito
  //   while(fgetc(fileLAS)!='\n');
  //  	// printf("[%u]  %.2lf, %.2lf, %.2lf\n",pointer[i].id,pointer[i].x,pointer[i].y,pointer[i].z);
  // }
  printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",xmin,xmax,ymin,ymax );


  double Width = roundl((xmax-xmin)*100)/100;
  double High = roundl((ymax-ymin)*100)/100;
  double Density = Npoints/(Width*High);
  printf("Ancho:  %.3lf\n",Width);
  printf("Alto:  %.3lf\n",High);
  printf("Densidad:  %.3lf\n",Density);


  unsigned short minNumPoints = 0.5*Density*Wsize*Wsize;
  printf("Numero minimo de puntos:   %u\n", minNumPoints);

  unsigned short Displace = Wsize*(1-Overlap);
  printf("Displacement   %u\n", Displace);
  unsigned short NumDisplace = Wsize/Displace-1;
  printf("N displacements   %u\n", NumDisplace);

  unsigned short Crow=(int)round(Width/Wsize), Ccol=(int)round(High/Wsize);
  printf("Celdas por columa: %d\n",Ccol);
  printf("Celdas por fila:   %d\n",Crow);

  unsigned short Ncells=Crow*Ccol;
  printf("N cells:    %u\n", Ncells);

  char command[255];
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



  double xxmin,xxmax,yymin,yymax,zzmin,zzvalue;
  // printf("Puntos por celda  :");
  int loopX, loopY=1, cellNum=0, cellPoints=0;
  yymin=ymin;
  yymax=ymin+Wsize;
  while(yymax<ymax){
      xxmin=xmin;
      xxmax=xmin+Wsize;
      loopX=1;
      while(xxmax<xmax){
          // sprintf(command,"las2las -i ../datos/sample24.las -o ../datos/window.las --minx %f --miny %f --maxx %f --maxy %f",xxmin,yymin,xxmax,yymax);
          sprintf(command,"las2las -i ../datos/sample24.las -o ../datos/window.las -e \"%f %f %f %f\"",xxmin,yymin,xxmax,yymax);
          if(system(command) < 0){
            perror("system");
            return 1;
          }
          reader = LASReader_Create("../datos/window.las");
          header = LASReader_GetHeader(reader);
          cellPoints = LASHeader_GetPointRecordsCount(header);
          if(cellPoints >= minNumPoints ){
              p = LASReader_GetNextPoint(reader);
              zzmin=LASHeader_GetMinZ(header);
              // printf("Minimo de la celda %d: %.2f\n",cellNum+1,zzmin);
              //Busco el mínimo; tengo en cuenta los problemas numericos de los double
              while(round(LASPoint_GetZ(p)*100)/100 > zzmin){
                // if(cellNum==54) printf("Z: %.4f, mayor? %d\n", LASPoint_GetZ(p),LASPoint_GetZ(p) > zzmin);
                p = LASReader_GetNextPoint(reader);
              }
              printf("[%d , %.2f]\n ",cellPoints,LASPoint_GetZ(p));
          }
          loopX++;
          xxmin=xmin+(loopX-1)*Wsize;
          xxmax=xmin+loopX*Wsize;
          cellNum++;

      }
      loopY++;
      yymin=ymin+(loopY-1)*Wsize;
      yymax=ymin+loopY*Wsize;
  }
  printf("\nCeldas      %d\n", cellNum);

  printf("\n\n/////////////////////////// END ///////////////////////////\n");

  // Compruebo que mínimos se repiten
  if(NumDisplace > 0){

  }

  // Aplico la malla para ver si hay zonas sin puntos
  if(Bsize > 0){
    printf("Aplico la malla\n");
  }
  // Libero memoria y cierro el fichero
  // free(cloud);
  // if(fclose(fileLAS)){
  //   printf("Cannot close the file\n");
  //   return -1;
  // }
  printf("Elapsed:     %.6f s\n",(double)(clock()-t_init) / CLOCKS_PER_SEC );
  return 0;
}
