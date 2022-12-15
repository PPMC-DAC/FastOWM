// gcc OWM.c -o func_OWM.o

#include <stdio.h>
#include <stdlib.h>
// #include <liblas/capi/las_version.h> //porque esta se llama ya desde liblas.h
// #include <liblas/capi/liblas.h>

int main(){

  typedef struct{
    double X, Y, Z;
  }Lpoint;

  // unsigned int discard;
  char discard;

  // LASpoint* puntero;

  FILE* ficheroLAS;
  if((ficheroLAS = fopen("../datos/sample54.xyz","r")) == NULL){
    printf("Unable to open file!\n");
		return -1;
  }

  // Recorre todo el fichero
  // while(!feof(ficheroLAS)) {
  int ch,lines=0;
  // int ch;
  while(!feof(ficheroLAS))
  {
    ch = fgetc(ficheroLAS);
    if(ch == '\n')
      lines++;
  }
  rewind(ficheroLAS);

  printf("%d lines\n",lines );


  Lpoint* puntero = malloc(lines* sizeof(Lpoint));
  for(int i=0; i<lines ; i++){
    //Obtengo los datos XYZ
    fscanf(ficheroLAS, "%lf%lf%lf" ,&puntero[i].X,&puntero[i].Y,&puntero[i].Z);
    //Los siguientes 6 campos no los necesito
    while(fgetc(ficheroLAS)!='\n');
   	printf("[%d]  %.2lf, %.2lf, %.2lf\n",i,puntero[i].X,puntero[i].Y,puntero[i].Z);
  }

  fclose(ficheroLAS);
  return 0;
}
