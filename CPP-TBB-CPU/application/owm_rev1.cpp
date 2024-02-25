// Difference wrt cesgaOWM.c 
// 1. Moves from C to C++
// 2. Octree replaced by Quadtree, that is now a class with its constructor
// 3. MinRadius (min size of a leaf bounding box) is now an imput parameter and argument of insertPoint instead of a constant

// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <time.h>
// #include <unistd.h>
// #include <string.h>
// #include <omp.h>

#include "../include/rev1_func.h"

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int main( int argc, char* argv[]){

    // Ficheros
    FILE* fileXYZ;
    FILE* fileMin;
#ifdef DEBUG
    FILE* fileDeb1, fileDeb2;
#endif
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

    unsigned int numLLPs =0;

    unsigned int addMin =0;

    double t_stage, t_func, t_quadtree;

    // Tamaño de la ventana deslizante
    unsigned short Wsize = 10;
    // Tamaño de la rejilla
    unsigned short Bsize = 20;
    // Solape de la ventana deslizante
    double Overlap = 0.8;
    // Numero de procesadores
    unsigned short num_procs = 4;

    //Control del bucle de ejecución
    unsigned int numRuns = 1;
    // chunk size used in the dynamic scheduling of the stage 1
    unsigned int chunk = 1;

    char inputTXT[128] = {"./data/INAER_2011_Alcoy.xyz"};
    // char inputLAS[50] = {"./data/INAER_2011_Alcoy.las"};
    char outputTXT[128] = {"./data/INAER_2011_Alcoy_salida.xyz"};
    // char outputLAS[50] = {"./data/INAER_2011_Alcoy_salida.las"};

    // Compruebo los argumentos de entrada
    if(argc>1) {
      strcpy(inputTXT,argv[1]);
      strcat(inputTXT,".xyz");

      strcpy(outputTXT,argv[1]);
      strcat(outputTXT,"_salida.xyz");
    }
    if(argc>2) Wsize = atoi(argv[2]);
    if(argc>3) Bsize = atoi(argv[3]);
    if(argc>4) Overlap = atof(argv[4]);
    if(argc>5) num_procs = atoi(argv[5]);
    if(argc>6) numRuns = atoi(argv[6]);
    // chunk size used in the dynamic scheduling of the stage 1
    if(argc>7) chunk = atoi(argv[7]);

    // qtree radius
    float minRadius = 0.1;

    double* resultados=new double[numRuns];

    omp_set_num_threads(num_procs);

    printf("Input.txt: %s ---> EX. CON %d CORES\n", inputTXT, num_procs);

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
    } else if(!strcmp(inputTXT,"./data/sample24.xyz")){
      Npoints = 7492;
      min.x   = 513748.12;
      max.x   = 513869.97;
      min.y   = 5403124.76;
      max.y   = 5403197.20;
    } else {// For files with header values (Npoints, min.x, max.x, min.y, max.y)
    //Read header values
      if(fscanf(fileXYZ, "%d\n%lf\n%lf\n%lf\n%lf\n",&Npoints, &min.x, &max.x, &min.y, &max.y) < 5){
        printf("Imposible to read header values\n");
        exit(-1);
        }
    }

    // Reservo memoria para la nube de puntos
    pointer = (Lpoint*)malloc(Npoints*sizeof(Lpoint));

    printf("VOLCANDO PUNTOS...\n");

    for(int i=0; i<Npoints ; i++){
      //Obtengo los datos id X Y Z
      pointer[i].id = i;
      if(fscanf(fileXYZ, "%lf %lf %lf",&pointer[i].x,&pointer[i].y,&pointer[i].z) < 3){
        printf("Imposible to obtain values\n");
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

    radius = getRadius(min, max, &maxRadius);
    center = getCenter(min, radius);
    printf("QTREE PARAMETERS:\n");
    printf("MaxRadius:  %.3f\n", maxRadius);
    printf("Center:     %.2f , %.2f\n", center.x,center.y);
    printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
    printf("CREANDO QTREE...\n");
    // qtreeIn = createQtree(center, maxRadius);
    Vector2D newCenter = {center.x,center.y};

    // Inserto los puntos en el Qtree
    printf("Inserting points; minRadius: %g\n", minRadius);

    t_quadtree=omp_get_wtime();
    Qtree qtreeIn = new Qtree_t( newCenter, maxRadius );
    for(int i = 0; i < Npoints; i++)
       insertPointF(&pointer[i], qtreeIn, minRadius);
    printf("Time elapsed at Quadtree construction:     %.6f s\n\n", omp_get_wtime()-t_quadtree);

    Width = round2d(max.x-min.x);
    High = round2d(max.y-min.y);
    // Densidad en puntos/m²
    Density = Npoints/(Width*High);
    printf("CLOUD PARAMETERS:\n");
    printf("Número de puntos      %d\n",Npoints);
    printf("Ancho:  %.2lf\n",Width);
    printf("Alto:  %.2lf\n",High);
    printf("Densidad:  %.3lf\n",Density);

    printf("\nTamaño de ventana     %u\n", Wsize);
    printf("Tamaño de rejilla     %u\nSolape                %.2f\n", Bsize,Overlap);

    // El numero minimo sera la mitad del numero de puntos medio por celda
    minNumPoints = 0.5*Density*Wsize*Wsize;
    printf("Minimo numero de puntos por celda:   %u\n", minNumPoints);

    Displace = round2d(Wsize*(1-Overlap));
    printf("Displacement   %.2f\n", Displace);

    // Stage 1 parameters
    printf("\nVENTANA:\n");

    if(Overlap > 0.0) {
     Crow=(int)(round((Width+2*Wsize*Overlap)/Displace))-1;
     Ccol=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
    } else {
     Crow=(int)floor(Width/Wsize)+1;
     Ccol=(int)floor(High/Wsize)+1;
    }
    printf("Celdas por columa: %d\n",Ccol);
    printf("Celdas por fila:   %d\n",Crow);
    Ncells = Crow*Ccol;

    // Stage 3 parameter
    printf("\nMALLA:\n");
    Crowg=(int)floor(Width/Bsize)+1;
    Ccolg=(int)floor(High/Bsize)+1;
    // Crow=(int)round(Width/Bsize);
    // Ccol=(int)round(High/Bsize);
    printf("Dimensiones de la malla %dx%d\n\n", Ccolg, Crowg);
    Ngrid = Crowg*Ccolg;
    // Voy a tener como mucho un mínimo por rejilla..
    minIDs = (int*)malloc(Ncells*sizeof(int));
    minGridIDs = (int*)malloc(Ngrid*sizeof(int));

    while(numRuns){

        t_func=omp_get_wtime();
        t_stage=omp_get_wtime();

        #ifdef TASKS
        countMin = stage1t(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min, chunk);
        #else
        countMin = stage1(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min, chunk);
        #endif

        printf("Time elapsed at STAGE 1:     %.6f s\n\n", omp_get_wtime()-t_stage);

        numLLPs=countMin;

        if(Overlap != 0){
            t_stage=omp_get_wtime();
            qsort(minIDs,countMin,sizeof(int),&cmpfunc);

#ifdef DEBUG
          if((fileDeb1 = fopen("sortedmins.txt","w")) == NULL){
            printf("Unable to create file!\n");
            return -1;
          }
          for(int i=0 ; i<countMin ; i++)
              fprintf(fileDeb1, "%d %.2f %.2f %.15f\n", minIDs[i], pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);
          fclose(fileDeb1);
#endif
            // Detect repeated ids and store them in minIDs. NumLLPs is the number of LLPs found and stored at the beggining of minIDs
                  numLLPs = stage2(countMin, minIDs);
#ifdef DEBUG
          if((fileDeb1 = fopen("LLPs.txt","w")) == NULL){
            printf("Unable to create file!\n");
            return -1;
          }
          for(int i=0 ; i<numLLPs ; i++)
              fprintf(fileDeb1, "%d %.2f %.2f %.15f\n", minIDs[i], pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);
          fclose(fileDeb1);
#endif            
            printf("Time elapsed at STAGE 2:     %.6f s\n\n",omp_get_wtime() - t_stage );
        }

        if(Bsize > 0){
            t_stage=omp_get_wtime();

            Qtree grid =  new Qtree_t( newCenter, maxRadius) ;
            for(int i = 0; i < numLLPs; i++)
               insertPointF(&pointer[minIDs[i]], grid, minRadius);
            addMin = stage3(Bsize, Crowg, Ccolg, minGridIDs, qtreeIn, grid, min);
            printf("Time elapsed at STAGE 3:     %.6f s\n\n",omp_get_wtime() - t_stage );
            deleteQtree(grid);
            delete(grid);
        }

        printf("TOTAL time elapsed:     %.6f s\n", resultados[--numRuns] = omp_get_wtime() - t_func);
        printf("Output ground seed-point cloud with %d points, %d fewer points than input cloud\n", numLLPs+addMin, Npoints - numLLPs - addMin );

        printf("/////////////////////////////////////////////////////////////////////\n");
    }

    printf("Time of each run:  ");
    printf("  %.4lf  ", resultados[0]);
    numRuns = atoi(argv[6]);
    if(numRuns > 0){
      for( int i=1 ; i<numRuns ; i++ ){
        printf("  %.4lf  ", resultados[i]);
        resultados[0] += resultados[i];
      }
      printf("\nAverage: %.4lf\n\n", resultados[0]/numRuns);
    } else 
      printf("\nAverage: %.4lf\n\n", resultados[0]);

#ifdef DEBUG
    // Fichero de salida
    printf("Creo el fichero %s ...\n", outputTXT);
    if((fileMin = fopen(outputTXT,"w")) == NULL){
      printf("Unable to create file!\n");
      return -1;
    }

    for(int i=0 ; i<numLLPs ; i++)
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
    // free(qtreeIn);
    // qtreeIn=NULL;
    free(minIDs);
    minIDs=NULL;
    free(minGridIDs);
    minGridIDs=NULL;
    free(pointer);
    pointer=NULL;

    deleteQtree(qtreeIn);
    delete(qtreeIn);


    return 0;
}
