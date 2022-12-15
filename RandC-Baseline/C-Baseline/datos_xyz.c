// gcc OWM.c -o func_OWM.o

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <liblas/capi/las_version.h> //porque esta se llama ya desde liblas.h
// #include <liblas/capi/liblas.h>

int main(int argc, char* argv[]){

  typedef struct{
    unsigned int id;
    double x;
    double y;
    double z;
  } Lpoint;

  // unsigned int discard;
  char discard;

  // LASpoint* puntero;

  char inputTXT[50] = {"../datos/INAER_2011_Alcoy.xyz"};
  // char inputLAS[50] = {"../datos/INAER_2011_Alcoy.las"};


  // Compruebo los argumentos de entrada
  if(argc>1) {
    strcpy(inputTXT,argv[1]);
  }

  FILE* fileLAS;
  if((fileLAS = fopen(inputTXT,"r")) == NULL){
    printf("Unable to open file!\n");
		return -1;
  }

  // Recorre todo el fichero
  // while(!feof(ficheroLAS)) {
  int ch,Npoints=0;
  // int ch;
  while(!feof(fileLAS))
  {
    ch = fgetc(fileLAS);
    if(ch == '\n')
      Npoints++;
  }
  rewind(fileLAS);

  printf("Numero de puntos: %d\n",Npoints );


  // Reservo memoria para todos los puntos y obtengo los datos junto los max y min
  double auxX, auxY, auxZ;
  Lpoint min, max;
  //Obtengo los datos id X Y Z
  // pointer[0].id = 0;
  fscanf(fileLAS, "%lf%lf%lf",&auxX,&auxY,&auxZ);
  // Los inicializo así para que no haya problemas con los mínimos
  min.x=auxX;
  max.x=auxX;
  min.y=auxY;
  max.y=auxY;
  //Los siguientes 6 campos no los necesito
  while(fgetc(fileLAS)!='\n');

  for(int i=1; i<Npoints-1 ; i++){
    //Obtengo los datos id X Y Z
    // pointer[i].id = i;
    fscanf(fileLAS, "%lf%lf%lf",&auxX,&auxY,&auxZ);
    if(min.x>auxX) min.x=auxX;
    if(max.x<auxX) max.x=auxX;
    if(min.y>auxY) min.y=auxY;
    if(max.y<auxY) max.y=auxY;
    //Los siguientes 6 campos no los necesito
    while(fgetc(fileLAS)!='\n');
   	// printf("[%u]  %.2lf, %.2lf, %.2lf\n",pointer[i].id,pointer[i].x,pointer[i].y,pointer[i].z);
  }
  printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );

  fclose(fileLAS);
  return 0;
}
