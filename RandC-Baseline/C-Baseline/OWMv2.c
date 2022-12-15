// gcc OWM.c -o func_OWM.o -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <liblas/capi/las_version.h> //porque esta se llama ya desde liblas.h
// #include <liblas/capi/liblas.h>

int main( int argc, char* argv[]){

  typedef struct{
    unsigned int id;
    double x;
    double y;
    double z;
  } Lpoint;

  // Tamaño de la ventana deslizante
  unsigned short Wsize = 12;
  // Tamaño de la rejilla
  unsigned char Bsize = 20;
  // Solape de la ventana deslizante
  float Overlap = 0;


  // Compruebo los argumentos de entrada
  if(argc>1) Wsize = atoi(argv[1]);
  if(argc>2) Bsize = atoi(argv[2]);
  if(argc>3) Overlap = atoi(argv[3]);

  printf("PARAMETERS:\nTamaño de ventana     %u\n", Wsize);
  printf("Tamaño de rejilla     %u\nSolape                %.2f\n", Bsize,Overlap);

  // Abro el fichero
  FILE* fileLAS;
  if((fileLAS = fopen("../datos/sample24.xyz","r")) == NULL){
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

  // Reservo memoria para todos los puntos y obtengo los datos junto los max y min
  double xmax=0, xmin=0, ymax=0, ymin=0;
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

  unsigned short minNumPoints = 0.5*Density*Wsize*Wsize;
  printf("Numero minimo de puntos:   %u\n", minNumPoints);

  unsigned short Displace = Wsize*(1-Overlap);
  printf("Displacement   %u\n", Displace);
  unsigned short NumDisplace = Wsize/Displace-1;
  printf("N displacements   %u\n", NumDisplace);

  unsigned short Ncells=round((Width/Wsize)*(High/Wsize));
  printf("N cells:    %u\n", Ncells);

  printf("/////////////////////////// LOOP ///////////////////////////\n");

  int j,i;
  for(j=-NumDisplace; j<=NumDisplace; j++ ){
    for(i=-NumDisplace; i<=NumDisplace; i++){
      printf("%d %d\n",j,i );
    }
  }


  printf("/////////////////////////// END ///////////////////////////\n");

  // Libero memoria y cierro el fichero
  free(pointer);
  if(fclose(fileLAS)){
    printf("Cannot close the file\n");
    return -1;
  }
  return 0;
}
