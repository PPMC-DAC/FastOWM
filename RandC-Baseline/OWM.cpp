// gcc OWM.c -o func_OWM.o

#include <liblas/liblas.hpp>

#include <fstream>  // std::ofstream
#include <algorithm> // std::copy
#include <exception> // std::exception

int main(){

  FILE *ficheroLAS;
  ficheroLAS = fopen("../datos/sample54.xyz","r");
  // Recorre todo el fichero
  // while(!feof(ficheroLAS)) fputc(fgetc(ficheroLAS), stdout);

  LASReaderH reader(ficheroLAS);

  fclose(ficheroLAS);
  return 0;
}
