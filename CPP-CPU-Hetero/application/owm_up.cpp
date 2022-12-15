
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <time.h>
// #include <unistd.h>
// #include <string.h>
// #include <omp.h>

#include "../include/envi_up.h"

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a < *(int*)b );
}

int main( int argc, char* argv[]){

    // Ficheros
    FILE* fileLAS;
    FILE* fileMin;

    // Octrees
    // Octree octreeIn = NULL;
    // Octree grid = NULL;

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

    double t_stage, t_func;

    // Tamaño de la ventana deslizante
    unsigned short Wsize = 12;
    // Tamaño de la rejilla
    unsigned short Bsize = 20;
    // Solape de la ventana deslizante
    double Overlap = 0.5;
    // Numero de procesadores
    unsigned short num_procs = 4;

    //Control del bucle de ejecución
    unsigned int bucle_entrada = 1;

    char inputTXT[128] = {"../datos/INAER_2011_Alcoy.xyz"};
    // char inputLAS[50] = {"../datos/INAER_2011_Alcoy.las"};
    char outputTXT[128] = {"../datos/INAER_2011_Alcoy_salida.xyz"};
    // char outputLAS[50] = {"../datos/INAER_2011_Alcoy_salida.las"};

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
    if(argc>6) bucle_entrada = atoi(argv[6]);

    double resultados[bucle_entrada];

    omp_set_num_threads(num_procs);

    printf("Input.txt: %s ---> EX. CON %d CORES\n", inputTXT, num_procs);

    // Abro el fichero
    if((fileLAS = fopen(inputTXT,"r")) == NULL){
      printf("Unable to open file!\n");
      exit(-1);
    }

    if( !strcmp(inputTXT,"../datos/INAER_2011_Alcoy.xyz") ){
      Npoints = 2772832;
      min.x   = 715244.96;
      max.x   = 716057.75;
      min.y   = 4286623.63;
      max.y   = 4287447.70;
      // min.z   = 0;
      // max.z   = 0; //No lo consulto nunca
    } else if( !strcmp(inputTXT,"../datos/INAER_2011_Alcoy_Core.xyz") ){
      Npoints = 20380212;
      min.x   = 714947.98;
      max.x   = 716361.06;
      min.y   = 4286501.93;
      max.y   = 4288406.23;
      // min.z   = 0;
      // max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/BABCOCK_2017_Arzua_3B.xyz")){
      Npoints = 40706503;
      min.x   = 568000.00;
      max.x   = 568999.99;
      min.y   = 4752320.00;
      max.y   = 4753319.99;
      // min.z   = 0;
      // max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/V21_group1_densified_point_cloud.xyz")){
      Npoints = 42384876;
      min.x   = 526964.093;
      max.x   = 527664.647;
      min.y   = 4742610.292;
      max.y   = 4743115.738;
      // min.z   = 0;
      // max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/V19_group1_densified_point_cloud.xyz")){
      Npoints = 48024480;
      min.x   = 526955.908;
      max.x   = 527686.445;
      min.y   = 4742586.025;
      max.y   = 4743124.373;
      // min.z   = 0;
      // max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/sample24.xyz")){
      Npoints = 7492;
      min.x   = 513748.12;
      max.x   = 513869.97;
      min.y   = 5403124.76;
      max.y   = 5403197.20;
      // min.z   = 0;
      // max.z   = 0; //No lo consulto nunca
    } else {
      printf("No header data!\n");
      exit(-1);
    }

    // Reservo memoria para la nube de puntos
    pointer = (Lpoint*)malloc(Npoints*sizeof(Lpoint));

    printf("VOLCANDO PUNTOS...\n");

    for(int i=0; i<Npoints ; i++){
      //Obtengo los datos id X Y Z
      pointer[i].id = i;
      if(fscanf(fileLAS, "%lf %lf %lf",&pointer[i].x,&pointer[i].y,&pointer[i].z) < 3){
        printf("Imposible to obtain values\n");
        exit(-1);
      }
      while(fgetc(fileLAS)!='\n');
    }

    //Ya no necesito mas el fichero
    if(fclose(fileLAS)){
      printf("Cannot close the file\n");
      exit(-1);
    }

    printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );

    radius = getRadius(min, max, &maxRadius);
    center = getCenter(min, radius);
    printf("OCTREE PARAMETERS:\n");
    printf("MaxRadius:  %.3f\n", maxRadius);
    printf("Center:     %.2f , %.2f\n", center.x,center.y);
    printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
    printf("CREANDO OCTREE...\n");
    // octreeIn = createOctree(center, maxRadius);
    Vector2D newCenter = {center.x,center.y};
    // Octree octreeIn = new Octree_t( newCenter, maxRadius );
    uOctree octreeIn( new nLeaf(newCenter, maxRadius) );

    // Inserto los puntos en el Octree
    for(int i = 0; i < Npoints; i++)
       insertPointF(&pointer[i], octreeIn, MIN_RADIUS);

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

    // sleep(2);
    while(bucle_entrada){

        // printf("\nN cells:    %u\n", Ncells);
        // Voy a tener como mucho un mínimo por ventana..
        // minIDs = malloc(Ncells*sizeof(int));
        t_func=omp_get_wtime();
        t_stage=omp_get_wtime();
        // printf("/////////////////////////// LOOP ///////////////////////////\n\n");


        // Me devuelve un mínimo por cada ventana no descartada y guarda el ID en minIDs
        // countMin = stage1(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, octreeIn.get(), min);
        countMin = stage1s(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, octreeIn.get(), min);

        printf("\nCeldas no descartadas:   %d\n", countMin);

        // printf("\n\n/////////////////////////// END ///////////////////////////\n");
        printf("Time elapsed at STAGE 1:     %.6f s\n\n", omp_get_wtime()-t_stage);

        // Para el caso de no hacer solpado; que searcher tenga un valor
        searcher=countMin;

        // Descarto mínimos si hay solape
        // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
        if(Overlap != 0){
            t_stage=omp_get_wtime();
            // printf("\nN minimos de entrada    %d\n", searcher);
            // printf("/////////////////////////// MIN SELECT ///////////////////////////\n\n");

            // Ordeno el array de IDs
            qsort(minIDs,countMin,sizeof(int),&cmpfunc);

            // for(int i=0 ; i<searcher ; i++)
            //   printf("%d %.2f\n", pointer[minIDs[i]].id,pointer[minIDs[i]].z);
            // Me quedo solo con los mínimos que se han repetido más de una vez
            searcher = stage2(countMin, minIDs);

            printf("\nNumero de minimos que me quedo: %d \n", searcher);

            // printf("\n\n/////////////////////////// END ///////////////////////////\n");
            printf("Time elapsed at STAGE 2:     %.6f s\n\n",omp_get_wtime() - t_stage );
        }


        // Aplico la malla para ver si hay zonas sin puntos
        if(Bsize > 0){
            // printf("N cells grid:    %u\n", Ngrid);
            // Voy a tener como mucho un mínimo por rejilla..
            // minGridIDs = malloc(Ngrid*sizeof(int));
            t_stage=omp_get_wtime();
            // printf("/////////////////////////// GRID ///////////////////////////\n\n\n");

            // Creo un nuevo octree con todos los mínimos; las mismas dimensiones que el grande
            // grid = createOctree(center, maxRadius);
            // Octree grid =  new Octree_t( newCenter, maxRadius) ;
            uOctree grid( new nLeaf(newCenter, maxRadius) );

            for(int i = 0; i < searcher; i++)
               insertPointF(&pointer[minIDs[i]], grid, MIN_RADIUS);

            //
            // addMin = stage3(Bsize, Crowg, Ccolg, minGridIDs, octreeIn.get(), grid.get(), min);
            addMin = stage3s(Bsize, Crowg, Ccolg, minGridIDs, octreeIn.get(), grid.get(), min);

            // printf("Minimos añadidos:        %d\n", addMin);
            // printf("\n\n/////////////////////////// END ///////////////////////////\n");
            printf("Time elapsed at STAGE 3:     %.6f s\n\n",omp_get_wtime() - t_stage );
            // Ya no necesito este octree
            deleteOctree(grid.get());
            // delete(grid);
            grid.reset(NULL);
        }

        printf("TOTAL time elapsed:     %.6f s\n", resultados[--bucle_entrada] = omp_get_wtime() - t_func);
        printf("Finalmente, el mapa va a tener %d puntos, %d puntos menos\n", searcher+addMin, Npoints - searcher+addMin );


        if(bucle_entrada){
          // free(minIDs);
          // minIDs=NULL;
          // free(minGridIDs);
          // minGridIDs=NULL;
          // Para parar un poco entre vuelta y vuelta
          sleep(20);
        }

        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");

    }

    printf("Ejecuciones:  ");
    for( int i=0 ; i<atoi(argv[6]) ; i++ ){
      printf("  %.4lf  ", resultados[i]);
      if(resultados[0] > resultados[i])
          resultados[0] = resultados[i];
    }
    printf("\nBEST: %.4lf\n", resultados[0]);

    // Fichero de salida
    printf("Creo el fichero %s ...\n", outputTXT);
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

    // Libero memoria
    // free(octreeIn);
    // octreeIn=NULL;
    free(minIDs);
    minIDs=NULL;
    free(minGridIDs);
    minGridIDs=NULL;
    free(pointer);
    pointer=NULL;

    deleteOctree(octreeIn.get());
    // delete(octreeIn);
    octreeIn.reset(NULL);


    return 0;
}
