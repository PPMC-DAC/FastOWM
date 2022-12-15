#include "file.h"



int readPoint(char *string, Lpoint *point){
  char * pch;
  //printf ("Splitting string \"%s\" into tokens:\n",string);
  pch = strtok (string," ,\t");
  if (pch == NULL) return -1;
  //printf("Retrieved %s\n", pch);
  //x
  point->x= atof(pch);
  //y
  pch = strtok (NULL, " ,\t");
  //printf("Retrieved %s\n", pch);
  if (pch == NULL) return -1;
  point->y= atof(pch);
  //z
  pch = strtok (NULL, " ,\t");
  //printf("Retrieved %s\n", pch);
  if (pch == NULL) return -1;
  point->z= atof(pch);
  //intensity
  pch = strtok (NULL, " ,\t");
  //printf("Retrieved %s\n", pch);
  if (pch == NULL) return -1;
  point->intensity= atoi(pch);
  return 0;
} 

int writePointsWithClass(Lpoint *points, int numPoints, char *fileName){
    long int i = 0, j = 0;
    FILE *file = NULL;
    printf("Writing output file %s...\n", fileName);
    file = fopen(fileName, "w");
    if(file == NULL){
        fprintf(file, "Error opening file to print: %s\n", fileName);
        return -1;
    }
    //fprintf(file, "X\t\t\tY\t\t\tZ\t\tI\t\tG\tT\tPID\tRID\n");

    for(i = 0; i < numPoints; i++)
    {    
        fprintf(file, "%.3lf\t%.3lf\t%.3lf\t%.3lf\t%u\n", points[i].x, points[i].y, points[i].z, (double) points[i].intensity, points[i].gId);
    }   

    fclose(file);
}

