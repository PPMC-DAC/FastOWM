// This version differs from optim1_qtree in the following ways:
// 1.- The function 'stage1' does not use the 'searchNeighbors2D' first to collect the points
//     inside the SW and later find the minimum if the number of points is large enough.
//     Instead, it uses findValidMin to find the minimum and the count of points inside the SW during the tree traversal.
// 2.- findValidMin does not check x,y coordinates of all points of a leaf. If the leaf's BBox is fully overlaped by the SW
//     x,y coordinates are not read. If the overlap is partial, x,y coordinate of all points have to be checked
// 3.- Stage1 in previous version used a induction variable to write the array of minimums (minIDs) so a critical section was
//     necessary for the parallel version. Now this induction is removed and minIDs can be written in parallel.

#include "../include/optim2_func.hpp"
#include <algorithm>

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
    
    int* minIDs = NULL;
    int* minGridIDs = NULL;

    unsigned int Npoints=0, Ncells, Ngrid;
    Vector2D center, radius, min, max;
    float maxRadius = 0.0;
    double Width, High, Density, Displace;
    unsigned short  minNumPoints, Crow, Ccol, Crowg, Ccolg;

    unsigned int countMin;

    unsigned int numLLPs = 0;

    unsigned int addMin = 0;

    double t_stage, t_func, t_qtree;

    // Tamaño de la ventana deslizante
    unsigned short Wsize = 10;
    // Tamaño de la rejilla
    unsigned short Bsize = 20;
    // Solape de la ventana deslizante
    double Overlap = 0.8;
    // Numero de procesadores
    unsigned short num_threads = 1;

    //Control del bucle de ejecución
    unsigned int numRuns = 1;

    char inputTXT[128] = {"./data/INAER_2011_Alcoy.xyz"};
    char outputTXT[128] = {"./data/INAER_2011_Alcoy_salida.xyz"};

    if(argc>1) {
      strcpy(inputTXT,argv[1]);
      strcat(inputTXT,".xyz");

      strcpy(outputTXT,argv[1]);
      strcat(outputTXT,"_salida.xyz");
    }
    if(argc>2) Wsize = atoi(argv[2]);
    if(argc>3) Bsize = atoi(argv[3]);
    if(argc>4) Overlap = atof(argv[4]);
    if(argc>5) num_threads = atoi(argv[5]);
    if(argc>6) numRuns = atoi(argv[6]);
    float minRadius = (argc>7)? atof(argv[7]) : 0.1;
    int level = (argc>8)? atoi(argv[8]) : 5;

    double* resultados=new double[numRuns];

    //omp_set_num_threads(num_threads);
    tbb::global_control c{tbb::global_control::max_allowed_parallelism, num_threads};

    printf("Input.txt: %s ---> EX. CON %d CORES\n", inputTXT, num_threads);

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
    Npoints++; // we insert at position 0 an artificial point to compare: {0.0, 0.0, std::numeric_limits<double>::max()}
    // Reservo memoria para la nube de puntos
    Lpoint* cloud = (Lpoint*)malloc(Npoints*sizeof(Lpoint));
    cloud[0]={0.0, 0.0, std::numeric_limits<double>::max()};

    printf("VOLCANDO PUNTOS...\n");
    for(int i=1; i<Npoints ; i++){
      //Obtengo los datos id X Y Z
      //cloud[i].id = i;
      if(fscanf(fileXYZ, "%lf %lf %lf",&cloud[i].x,&cloud[i].y,&cloud[i].z) < 3){
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

    printf("Inserting points with minRadius: %g\n", minRadius);
    t_qtree=omp_get_wtime();
    // Qtree qtreeIn = new Qtree_t( center, maxRadius );
    // for(int i = 1; i <= Npoints; i++)
    //    insertPoint(i, cloud, qtreeIn, minRadius);
    
    Qtree qtreeIn = parallel_qtree( level, center, maxRadius, cloud, Npoints, minRadius );

    printf("Time elapsed at Quadtree construction:     %.6f s\n\n", omp_get_wtime()-t_qtree);

    Width = round2d(max.x-min.x);
    High = round2d(max.y-min.y);
    Density = (Npoints-1)/(Width*High);
    printf("CLOUD PARAMETERS:\n");
    printf("Número de puntos      %d\n",Npoints-1);
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
    printf("Dimensiones de la malla %dx%d\n\n", Ccolg, Crowg);
    Ngrid = Crowg*Ccolg;
    minIDs = (int*)malloc(Ncells*sizeof(int));
    minGridIDs = (int*)malloc(Ngrid*sizeof(int));

    int numRuns_saved=numRuns;
    while(numRuns){

        t_func=omp_get_wtime();
        t_stage=omp_get_wtime();


        // The array minIDs that store the valid min of each cell/SW is initialized with -1
        // The cells/SWs without a valid minimum will keep the -1
        //std::fill(minIDs, minIDs+Ncells, -1);
        stage1(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, cloud, qtreeIn, min);

        printf("Time elapsed at STAGE 1:     %.6f s\n\n", omp_get_wtime()-t_stage);

        if(Overlap != 0){
            t_stage=omp_get_wtime();
            std::sort(oneapi::dpl::execution::par_unseq,minIDs,minIDs+Ncells); //std::execution::par, std::execution::par_unseq,

#ifdef DEBUG
          if((fileDeb1 = fopen("sortedmins.txt","w")) == NULL){
            printf("Unable to create file!\n");
            return -1;
          }
          for(int i=0 ; i<Ncells ; i++)
              fprintf(fileDeb1, "%d %.2f %.2f %.15f\n", minIDs[i], cloud[minIDs[i]].x, cloud[minIDs[i]].y,cloud[minIDs[i]].z);
          fclose(fileDeb1);
#endif
            // Detect repeated ids and store them in minIDs. NumLLPs is the number of LLPs found and stored at the beggining of minIDs
                  numLLPs = stage2(Ncells, minIDs);
#ifdef DEBUG
          if((fileDeb1 = fopen("LLPs.txt","w")) == NULL){
            printf("Unable to create file!\n");
            return -1;
          }
          for(int i=0 ; i<numLLPs ; i++)
              fprintf(fileDeb1, "%d %.2f %.2f %.15f\n", minIDs[i], cloud[minIDs[i]].x, cloud[minIDs[i]].y,cloud[minIDs[i]].z);
          fclose(fileDeb1);
#endif

            printf("Time elapsed at STAGE 2:     %.6f s\n\n",omp_get_wtime() - t_stage );
        }
        else{
          countMin = std::count_if(minIDs, minIDs+Ncells, [](int i){return i>=0;} );
          printf("Number of found minima:   %d\n\n", countMin);
          numLLPs=countMin;
        }

        printf("Number of found LLPs: %d \n\n", numLLPs);

        if(Bsize > 0){
            t_stage=omp_get_wtime();
            Qtree grid =  new Qtree_t( center, maxRadius) ;

            for(int i = 0; i < numLLPs; i++)
               insertPoint(minIDs[i], cloud, grid, minRadius);

            addMin = stage3(Bsize, Crowg, Ccolg, minGridIDs, cloud, qtreeIn, grid, min);
            printf("Time elapsed at STAGE 3:     %.6f s\n\n",omp_get_wtime() - t_stage );
            printf("Number of points added at stage 3: %d \n\n", addMin);
            deleteQtree(grid);
            delete(grid);
        }

        printf("TOTAL time elapsed:     %.6f s\n", resultados[--numRuns] = omp_get_wtime() - t_func);
        printf("Output ground seed-point cloud with %d points, %d fewer points than input cloud\n", numLLPs+addMin, Npoints-1 - numLLPs+addMin );

        printf("/////////////////////////////////////////////////////////////////////\n");
    }

    printf("Time of each run:  ");
    printf("  %.4lf  ", resultados[0]);
    numRuns = numRuns_saved;
    for( int i=1 ; i<numRuns ; i++ ){
        printf("  %.4lf  ", resultados[i]);
        resultados[0] += resultados[i];
    }
    printf("\nAverage: %.4lf\n\n", resultados[0]/numRuns);
    
#ifdef DEBUG
    // Fichero de salida
    printf("Creo el fichero %s ...\n", outputTXT);
    if((fileMin = fopen(outputTXT,"w")) == NULL){
      printf("Unable to create file!\n");
      return -1;
    }

    for(int i=0 ; i<numLLPs ; i++)
      fprintf(fileMin, "%.2f %.2f %.2f\n", cloud[minIDs[i]].x, cloud[minIDs[i]].y,cloud[minIDs[i]].z);

    for(int i=0 ; i<addMin ; i++)
      fprintf(fileMin, "%.2f %.2f %.2f\n", cloud[minGridIDs[i]].x, cloud[minGridIDs[i]].y,cloud[minGridIDs[i]].z);

    //Ya no necesito mas el fichero
    if(fclose(fileMin)){
      printf("Cannot close the file\n");
      return -1;
    }
#endif
    // Libero memoria
    free(minIDs);
    minIDs=NULL;
    free(minGridIDs);
    minGridIDs=NULL;
    free(cloud);
    cloud=NULL;

    deleteQtree(qtreeIn);
    delete(qtreeIn);


    return 0;
}