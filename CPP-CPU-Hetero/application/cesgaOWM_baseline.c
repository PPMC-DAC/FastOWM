
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

#include "../include/environment.h"

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b ); // Portable ascending order
}

int main( int argc, char* argv[]){

    // Ficheros
    FILE* fileXYZ;
    FILE* fileMin;

#ifdef DEBUG
    FILE* fileDeb1, fileDeb2;
#endif

    // Octrees
    Octree octreeIn = NULL;
    Octree grid = NULL;

    //Listas
    Lpoint* pointer = NULL;
    int* minIDs = NULL;
    int* minGridIDs = NULL;

    unsigned int Npoints=0, Ncells, Ngrid;
    Vector3D center, radius, min, max;
    float maxRadius = 0.0;
    double Width, High, Density, Displace;
    unsigned short  minNumPoints, Crow, Ccol, Crowg, Ccolg;

    unsigned int countMin;

    unsigned int searcher =0;

    unsigned int addMin =0;

    double t_stage, t_func, t_octree;

    // Tamaño de la ventana deslizante
    unsigned short Wsize = 10;
    // Tamaño de la rejilla
    unsigned short Bsize = 20;
    // Solape de la ventana deslizante
    double Overlap = 0.8;
    // Numero de procesadores
    unsigned short num_procs = 1;

    //Control del bucle de ejecución
    unsigned int numRuns = 1;

    char inputTXT[128] = {"./data/INAER_2011_Alcoy.xyz"};
    // char inputLAS[50] = {"./data/INAER_2011_Alcoy.las"};
    char outputTXT[128] = {"./data/INAER_2011_Alcoy_output.xyz"};
    // char outputLAS[50] = {"./data/INAER_2011_Alcoy_output.las"};

    // Compruebo los argumentos de entrada
    if(argc>1) {
      strcpy(inputTXT,argv[1]);
      strcat(inputTXT,".xyz");

      strcpy(outputTXT,argv[1]);
      strcat(outputTXT,"_output.xyz");
    }
    if(argc>2) Wsize = atoi(argv[2]);
    if(argc>3) Bsize = atoi(argv[3]);
    if(argc>4) Overlap = atof(argv[4]);
    if(argc>5) num_procs = atoi(argv[5]);
    if(argc>6) numRuns = atoi(argv[6]);

    double resultados[numRuns];

#ifdef PARALLEL
    omp_set_num_threads(num_procs);
    printf("Input.txt: %s ---> EX. with %d CORES\n", inputTXT, num_procs);
#else
    printf("Input.txt: %s ---> SEQUENTIAL\n", inputTXT);
#endif

    // Abro el fichero
    if((fileXYZ = fopen(inputTXT,"r")) == NULL){
      printf("Unable to open file!\n");
      exit(-1);
    }

    if( !strcmp(inputTXT,"./data/INAER_2011_Alcoy.xyz") ){
      Npoints = 2772832;
      min.x   = 715244.96;
      max.x   = 716057.75;
      min.y   = 4286623.63;
      max.y   = 4287447.70;
    } else if( !strcmp(inputTXT,"./data/INAER_2011_Alcoy_Core.xyz") ){
      Npoints = 20380212;
      min.x   = 714947.98;
      max.x   = 716361.06;
      min.y   = 4286501.93;
      max.y   = 4288406.23;
    } else if(!strcmp(inputTXT,"./data/BABCOCK_2017_Arzua_3B.xyz")){
      Npoints = 40706503;
      min.x   = 568000.00;
      max.x   = 568999.99;
      min.y   = 4752320.00;
      max.y   = 4753319.99;
    } else if(!strcmp(inputTXT,"./data/V21_group1_densified_point_cloud.xyz")){
      Npoints = 42384876;
      min.x   = 526964.093;
      max.x   = 527664.647;
      min.y   = 4742610.292;
      max.y   = 4743115.738;
    } else if(!strcmp(inputTXT,"./data/V19_group1_densified_point_cloud.xyz")){
      Npoints = 48024480;
      min.x   = 526955.908;
      max.x   = 527686.445;
      min.y   = 4742586.025;
      max.y   = 4743124.373;
    } else { // For files with header values (Npoints, min.x, max.x, min.y, max.y)
    //Read header values
      if(fscanf(fileXYZ, "%d\n%lf\n%lf\n%lf\n%lf\n",&Npoints, &min.x, &max.x, &min.y, &max.y) < 5){
        printf("Imposible to read header values\n");
        exit(-1);
        }
    }


    // Reservo memoria para la nube de puntos
    pointer = (Lpoint*)malloc(Npoints*sizeof(Lpoint));

    printf("Reading points...\n");

    for(int i=0; i<Npoints ; i++){
      //Obtengo los datos id X Y Z
      pointer[i].id = i;
      if(fscanf(fileXYZ, "%lf %lf %lf",&pointer[i].x,&pointer[i].y,&pointer[i].z) < 3){
        printf("Imposible to read values\n");
        exit(-1);
      }
      while(fgetc(fileXYZ)!='\n');
    }

    //Ya no necesito mas el fichero
    if(fclose(fileXYZ)){
      printf("Cannot close the file\n");
      exit(-1);
    }

    printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );
    // maxRadius identifies in which axis (x, y, z) the BBox has a maximum radius
    // In other words, it is the max(radius.x, radius.y, radius.z)
    radius = getRadius(min, max, &maxRadius);
    center = getCenter(min, radius);
    printf("OCTREE PARAMETERS:\n");
    printf("MaxRadius:  %.3f\n", maxRadius);
    printf("Center:     %.2f , %.2f\n", center.x,center.y);
    printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
    printf("BUILDING OCTREE...\n");

    t_octree=omp_get_wtime();
    octreeIn = createOctree(center, maxRadius);

    // Inserto los puntos en el Octree
    for(int i = 0; i < Npoints; i++)
       insertPoint(&pointer[i], octreeIn);

    printf("Time elapsed at Octree construction:     %.6f s\n\n", omp_get_wtime()-t_octree);
    Width = round2d(max.x-min.x);
    High = round2d(max.y-min.y);
    // Densidad en puntos/m²
    Density = Npoints/(Width*High);
    printf("CLOUD PARAMETERS:\n");
    printf("Number of points (cloud size): %d\n",Npoints);
    printf("Width:  %.2lf\n",Width);
    printf("Height:  %.2lf\n",High);
    printf("Density:  %.3lf\n",Density);

    printf("\nSize of sliding window (SW) %u\n", Wsize);
    printf("Grid size     %u\nOverlap                %.2f\n", Bsize,Overlap);

    // A minimum of a slide window, SW, is valid if there are at least minNumPoints inside the SW
    minNumPoints = 0.5*Density*Wsize*Wsize;
    printf("Minium number of points in the SW to be considered:   %u\n", minNumPoints);

    Displace = round2d(Wsize*(1-Overlap));
    printf("SW x and y Displacement   %.2f\n", Displace);

    // Stage 1 parameters
    printf("\nSliding Window (SW) parameters:\n");

    if(Overlap > 0.0) {
     Crow=(int)(round((Width+2*Wsize*Overlap)/Displace))-1;
     Ccol=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
    } else {
     Crow=(int)floor(Width/Wsize)+1;
     Ccol=(int)floor(High/Wsize)+1;
    }
    printf("Number of SWs per column: %d\n",Ccol);
    printf("Number of SW per row:   %d\n",Crow);
    Ncells = Crow*Ccol;
    printf("Total number of OWM steps (Crow x CCol):   %d\n",Ncells);

    // Stage 3 parameter
    printf("\nGrid (for stage 3) parameters:\n");
    Crowg=(int)floor(Width/Bsize)+1;
    Ccolg=(int)floor(High/Bsize)+1;

    printf("Grid dimesions in %dx%d boxes: %dx%d\n\n", Bsize, Bsize, Ccolg, Crowg);
    Ngrid = Crowg*Ccolg;

    // Voy a tener como mucho un mínimo por ventana..
    minIDs = calloc(Ncells,sizeof(int));
    // Voy a tener como mucho un mínimo por rejilla..
    minGridIDs = malloc(Ngrid*sizeof(int));

    while(numRuns){

        t_func=omp_get_wtime();
        t_stage=omp_get_wtime();

        // Me devuelve un mínimo por cada ventana no descartada y guarda el ID en minIDs
        countMin = stage1(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, octreeIn, min);

        printf("Time elapsed at STAGE 1:     %.6f s\n", omp_get_wtime()-t_stage);
        printf("Number of found minima:   %d\n\n", countMin);

        // Descarto mínimos si hay solape
        // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
        if(Overlap > 0.0){
            t_stage=omp_get_wtime();

            // Ordeno el array de IDs
            //qsort(minIDs,Ncells,sizeof(int),&cmpfunc); //This is not correct (it consider point with id=0 as one of the minima)
            qsort(minIDs,countMin,sizeof(int),&cmpfunc);
#ifdef DEBUG
          if((fileDeb1 = fopen("sortedmins.txt","w")) == NULL){
            printf("Unable to create file!\n");
            return -1;
          }
          for(int i=0 ; i<countMin ; i++)
              fprintf(fileDeb1, "%d %.2f %.2f %.2f\n", minIDs[i], pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);
          fclose(fileDeb1);
#endif
            // Me quedo solo con los mínimos que se han repetido más de una vez
            //searcher = stage2(Ncells, minIDs);//This is not correct (it consider point with id=0 as one of the minima)
            searcher = stage2(countMin, minIDs);
#ifdef DEBUG
          if((fileDeb1 = fopen("LLPs.txt","w")) == NULL){
            printf("Unable to create file!\n");
            return -1;
          }
          for(int i=0 ; i<searcher ; i++)
              fprintf(fileDeb1, "%d %.2f %.2f %.2f\n", minIDs[i], pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);
          fclose(fileDeb1);
#endif

            printf("Time elapsed at STAGE 2:     %.6f s\n",omp_get_wtime() - t_stage );
            printf("Number of found LLPs: %d \n\n", searcher);

        }
        else
            searcher = countMin;// Para el caso de no hacer solapado; que searcher tenga un valor

        // Aplico la malla para ver si hay zonas sin puntos
        if(Bsize > 0){

            t_stage=omp_get_wtime();

            // Creo un nuevo octree con todos los mínimos; las mismas dimensiones que el grande
            grid = createOctree(center, maxRadius);
            for(int i = 0; i < searcher; i++)
               insertPoint(&pointer[minIDs[i]], grid);
            addMin = stage3(Bsize, Crowg, Ccolg, minGridIDs, octreeIn, grid, min);

            printf("Time elapsed at STAGE 3:     %.6f s\n",omp_get_wtime() - t_stage );
            printf("Number of points added at stage 3: %d \n\n", addMin);

            // Ya no necesito este octree
            free(grid);
            grid = NULL;
        }

        printf("TOTAL time elapsed:     %.6f s\n", resultados[--numRuns] = omp_get_wtime() - t_func);
        printf("Output ground seed-point cloud with %d points, %d fewer points than input cloud\n", searcher+addMin, Npoints - searcher+addMin );


        // if(numRuns){
        //   sleep(5);
        // }

        printf("/////////////////////////////////////////////////////////////////////\n");
    }

    printf("Time of each run:  ");
    printf("  %.4lf  ", resultados[0]);
    numRuns = atoi(argv[6]);
    if(numRuns > 1){
      for( int i=1 ; i<numRuns ; i++ ){
        printf("  %.4lf  ", resultados[i]);
        resultados[0] += resultados[i];
      }
      printf("\nAverage: %.4lf\n\n", resultados[0]/numRuns);
    } else 
      printf("\nAverage: %.4lf\n\n", resultados[0]);

#ifdef DEBUG
    // Output file
    printf("Open file %s ...\n", outputTXT);
    if((fileMin = fopen(outputTXT,"w")) == NULL){
      printf("Unable to create file!\n");
      return -1;
    }

    for(int i=0 ; i<searcher ; i++)
      fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);

    for(int i=0 ; i<addMin ; i++)
      fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minGridIDs[i]].x, pointer[minGridIDs[i]].y,pointer[minGridIDs[i]].z);

    //Ya no necesito mas el fichero
    if(fclose(fileMin)){
      printf("Cannot close the file\n");
      return -1;
    }
#endif
    // Libero memoria
    free(octreeIn);
    octreeIn=NULL;
    free(minIDs);
    minIDs=NULL;
    free(minGridIDs);
    minGridIDs=NULL;
    free(pointer);
    pointer=NULL;


    return 0;
}
