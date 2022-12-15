// gcc OWM.c -o func_OWM.o -lm -llas_c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <liblas/capi/las_version.h> //porque esta se llama ya desde liblas.h
#include <liblas/capi/liblas.h>

int main( int argc, char* argv[]){

  typedef struct{
    unsigned int id;
    double x;
    double y;
    double z;
  } Lpoint;

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
  FILE* fileLAS;
  if((fileLAS = fopen("../datos/sample24.xyz","r")) == NULL){
  // if((fileLAS = fopen("../datos/INAER_2011_Alcoy.xyz","r")) == NULL){
    printf("Unable to open file!\n");
		return -1;
  }

  // Recorro todo el file para contar las líneas
  unsigned int Npoints=0;
  while(!feof(fileLAS))
    if(fgetc(fileLAS) == '\n')
      Npoints++;
  printf("Número de puntos      %d\n",Npoints );

  // Vuelvo al principio del fichero
  rewind(fileLAS);

  // Reservo memoria para todos los puntos y obtengo los datos junto con los max y min
  double xmax=0, xmin=0, ymax=0, ymin=0;
  // Lpoint xmin, xmax, ymax, ymin;
  Lpoint* pointer = malloc(Npoints* sizeof(Lpoint));

  //Obtengo los datos id X Y Z
  pointer[0].id = 0;
  fscanf(fileLAS, "%lf%lf%lf",&pointer[0].x,&pointer[0].y,&pointer[0].z);
  // Los inicializo así para que no haya problemas con los mínimos
  xmin=pointer[0].x;
  xmax=pointer[0].x;
  ymin=pointer[0].y;
  ymax=pointer[0].y;
  //Los siguientes 6 campos no los necesito
  while(fgetc(fileLAS)!='\n');

  for(int i=1; i<Npoints-1 ; i++){
    //Obtengo los datos id X Y Z
    pointer[i].id = i;
    fscanf(fileLAS, "%lf%lf%lf",&pointer[i].x,&pointer[i].y,&pointer[i].z);
    if(xmin>pointer[i].x) xmin=pointer[i].x;
    if(xmax<pointer[i].x) xmax=pointer[i].x;
    if(ymin>pointer[i].y) ymin=pointer[i].y;
    if(ymax<pointer[i].y) ymax=pointer[i].y;
    //Los siguientes 6 campos no los necesito
    while(fgetc(fileLAS)!='\n');
   	// printf("[%u]  %.2lf, %.2lf, %.2lf\n",pointer[i].id,pointer[i].x,pointer[i].y,pointer[i].z);
  }
  printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",xmin,xmax,ymin,ymax );


  double Width = roundl((xmax-xmin)*100)/100;
  double High = roundl((ymax-ymin)*100)/100;
  double Density = Npoints/(Width*High);
  printf("Ancho:  %.3lf\n",Width);
  printf("Alto:  %.3lf\n",High);
  printf("Densidad:  %.3lf\n",Density);

  // Busco el numero de filas y columnas
  // Una diferencia grande para no tener que comprobarla siempre
  // double dif=pointer[0].x-pointer[(int)(Width/2)].x;
  // unsigned int row=0, init=0, index=(int)(Width/2+2);
  // unsigned int RowCol[100]={0};
  // // eps es para no tener que recorrer toda la fila; un margen
  // // TODO: Puede ser ajustable
  // unsigned int eps =1; //(int)(Width/2);
  // while(index < Npoints){
  //   printf("%u\n",index);
  //   while(dif<(pointer[init].x-pointer[index].x)){
  //     index++;
  //     //Compruebo que no me salgo y si lo voy a hacer, salgo del bucle
  //     if(index==Npoints) init=index;
  //     // printf("%u\n",index);
  //   }
  //   // Si llego aquí es porque tengo que cambiar de fila
  //   // Las filas son las posiciones del vector y las columnas el valor
  //   RowCol[row]=(index-init)-1;
  //   row++;
  //   init=index;
  // }
  // for(int i=0; i<100; i++) printf("%u  ",RowCol[i] );
  // printf("\n");

  unsigned short minNumPoints = 0.5*Density*Wsize*Wsize;
  printf("Numero minimo de puntos:   %u\n", minNumPoints);

  unsigned short Displace = Wsize*(1-Overlap);
  printf("Displacement   %u\n", Displace);
  unsigned short NumDisplace = Wsize/Displace-1;
  printf("N displacements   %u\n", NumDisplace);

  unsigned short Ncells=round((Width/Wsize)*(High/Wsize));
  printf("N cells:    %u\n", Ncells);

  printf("/////////////////////////// LOOP ///////////////////////////\n\n\n");

  int j,i;
  Lcell myCell = {0.0,0.0,0.0,0.0,0};
  for(j=-NumDisplace; j<=NumDisplace; j++ ){
    for(i=-NumDisplace; i<=NumDisplace; i++){
      printf("Point [%d, %d]\n",j,i );
      // myCell.xmin = pointer[Wsize];
      // myCell.xmax = pointer[0];
    }
  }


  printf("\n\n/////////////////////////// END ///////////////////////////\n");

  // Libero memoria y cierro el fichero
  free(pointer);
  if(fclose(fileLAS)){
    printf("Cannot close the file\n");
    return -1;
  }
  printf("Elapsed:     %.3f s\n",(float)(clock()-t_init) / CLOCKS_PER_SEC );
  return 0;
}
