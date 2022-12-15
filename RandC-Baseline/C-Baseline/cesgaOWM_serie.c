
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

#include "include/environment.h"

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int main( int argc, char* argv[]){

    // Ficheros
    FILE* fileLAS;
    FILE* fileMin;

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

    unsigned int searcher;

    unsigned int addMin;

    double t_stage, t_func;

    // Tamaño de la ventana deslizante
    unsigned short Wsize = 12;
    // Tamaño de la rejilla
    unsigned short Bsize = 20;
    // Solape de la ventana deslizante
    double Overlap = 0.8;
    // Numero de procesadores
    unsigned short num_procs = 4;

    //Control del bucle de ejecución
    unsigned int bucle_entrada = 1;

    char inputTXT[50] = {"../datos/INAER_2011_Alcoy.xyz"};
    // char inputLAS[50] = {"../datos/INAER_2011_Alcoy.las"};
    char outputTXT[50] = {"../datos/INAER_2011_Alcoy_salida.xyz"};
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

    // omp_set_num_threads(num_procs);

    printf("Input.txt: %s ---> EX. CON %d CORES\n", inputTXT, num_procs);

    // Abro el fichero
    if((fileLAS = fopen(inputTXT,"r")) == NULL){
      printf("Unable to open file!\n");
      return -1;
    }

    // // Numero de puntos para hacer la reserva de memoria justa
    // while(!feof(fileLAS))
    //   if(fgetc(fileLAS) == '\n')
    //     Npoints++;
    // // Vuelvo al principio del fichero
    // rewind(fileLAS);

    if( !strcmp(inputTXT,"../datos/INAER_2011_Alcoy.xyz") ){
      Npoints = 2772832;
      min.x   = 715244.96;
      max.x   = 716057.75;
      min.y   = 4286623.63;
      max.y   = 4287447.70;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if( !strcmp(inputTXT,"../datos/INAER_2011_Alcoy_Core.xyz") ){
      Npoints = 20380212;
      min.x   = 714947.98;
      max.x   = 716361.06;
      min.y   = 4286501.93;
      max.y   = 4288406.23;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/BABCOCK_2017_Arzua_3B.xyz")){
      Npoints = 40706503;
      min.x   = 568000.00;
      max.x   = 568999.99;
      min.y   = 4752320.00;
      max.y   = 4753319.99;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/INSITU_2018_Jenaro.xyz")){
      Npoints = 42978528;
      min.x   = 536067.12;
      max.x   = 536213.34;
      min.y   = 4746607.46;
      max.y   = 4746950.19;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/V21_group1_densified_point_cloud.xyz")){
      Npoints = 42384876;
      min.x   = 526964.09;
      max.x   = 527664.65;
      min.y   = 4742610.29;
      max.y   = 4743115.74;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/MMCoruna_020_Aj_Lp.xyz")){
      Npoints = 53467933;
      min.x   = 547248.67;
      max.x   = 547590.67;
      min.y   = 4801345.51;
      max.y   = 4801855.17;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/Begonte.xyz")){
      Npoints = 31797424;
      min.x   = 605523.26;
      max.x   = 606023.26;
      min.y   = 4782855.38;
      max.y   = 4783355.37;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/Guitiriz_foto.xyz")){
      Npoints = 27445172;
      min.x   = 591000.00;
      max.x   = 591500.00;
      min.y   = 4784500.01;
      max.y   = 4785000.00;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/Guitiriz.xyz")){
      Npoints = 10558049;
      min.x   = 591000.00;
      max.x   = 591500.00;
      min.y   = 4784500.01;
      max.y   = 4785000.00;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/sample24.xyz")){
      Npoints = 7492;
      min.x   = 513748.12;
      max.x   = 513869.97;
      min.y   = 5403124.76;
      max.y   = 5403197.20;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/VaihingenS03.xyz")){
      Npoints = 3776182;
      min.x   = 496400.74;
      max.x   = 497852.64;
      min.y   = 5418996.49;
      max.y   = 5419503.44;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/VaihingenS05.xyz")){
      Npoints = 3582656;
      min.x   = 496398.92;
      max.x   = 497850.53;
      min.y   = 5419290.41;
      max.y   = 5419793.69;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/VaihingenS07.xyz")){
      Npoints = 3728882;
      min.x   = 496400.06;
      max.x   = 497850.31;
      min.y   = 5418705.66;
      max.y   = 5419211.81;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/VaihingenS09.xyz")){
      Npoints = 3675745;
      min.x   = 496400.45;
      max.x   = 497851.11;
      min.y   = 5419580.33;
      max.y   = 5420081.71;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    } else if(!strcmp(inputTXT,"../datos/VaihingenS10.xyz")){
      Npoints = 3801516;
      min.x   = 496400.00;
      max.x   = 497850.28;
      min.y   = 5419886.80;
      max.y   = 5420386.40;
      min.z   = 0;
      max.z   = 0; //No lo consulto nunca
    }

    // Reservo memoria para la nube de puntos
    pointer = malloc(Npoints*sizeof(Lpoint));

    printf("VOLCANDO PUNTOS...\n");
    //Obtengo los datos id X Y Z
    // pointer[0].id = 0;
    // fscanf(fileLAS, "%lf%lf%lf",&pointer[0].x,&pointer[0].y,&pointer[0].z);
    // // Los inicializo así para que no haya problemas con los mínimos
    // min.x=pointer[0].x;
    // min.y=pointer[0].y;
    // //Los siguientes 6 campos no los necesito
    // while(fgetc(fileLAS)!='\n');

    for(int i=0; i<Npoints ; i++){
      //Obtengo los datos id X Y Z
      pointer[i].id = i;
      if(fscanf(fileLAS, "%lf%lf%lf",&pointer[i].x,&pointer[i].y,&pointer[i].z) < 3){
        printf("Imposible to obtain values\n");
        return -1;
      }
      // if(min.x>pointer[i].x) min.x=pointer[i].x;
      // if(max.x<pointer[i].x) max.x=pointer[i].x;
      // if(min.y>pointer[i].y) min.y=pointer[i].y;
      // if(max.y<pointer[i].y) max.y=pointer[i].y;
      //Los siguientes 6 campos no los necesito
      while(fgetc(fileLAS)!='\n');
     	// printf("[%u]  %.2lf, %.2lf, %.2lf\n",pointer[i].id,pointer[i].x,pointer[i].y,pointer[i].z);
    }

    //Ya no necesito mas el fichero
    if(fclose(fileLAS)){
      printf("Cannot close the file\n");
      return -1;
    }

    printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );

    radius = getRadius(min, max, &maxRadius);
    center = getCenter(min, radius);
    printf("OCTREE PARAMETERS:\n");
    printf("MaxRadius:  %.3f\n", maxRadius);
    printf("Center:     %.2f , %.2f\n", center.x,center.y);
    printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
    printf("CREANDO OCTREE...\n");
    octreeIn = createOctree(center, maxRadius);

    // Inserto los puntos en el Octree
    for(int i = 0; i < Npoints; i++)
       insertPoint(&pointer[i], octreeIn);


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

    printf("\nVENTANA:\n");
    // Crow=(int)round(Width/Wsize);
    // Ccol=(int)round(High/Wsize);
    Crow=(int)floor(Width/Wsize)+1;
    Ccol=(int)floor(High/Wsize)+1;
    if(Overlap!=0) Crow=(int)(round((Width+2*Wsize*Overlap)/Displace))-1, Ccol=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
    printf("Celdas por columa: %d\n",Ccol);
    printf("Celdas por fila:   %d\n",Crow);
    Ncells = Crow*Ccol;

    // Voy a tener como mucho un mínimo por ventana..
    minIDs = malloc(Ncells*sizeof(int));

    // PARAMETROS DE LA MALLA
    printf("\nMALLA:\n");
    Crowg=(int)floor(Width/Bsize)+1;
    Ccolg=(int)floor(High/Bsize)+1;
    // Crow=(int)round(Width/Bsize);
    // Ccol=(int)round(High/Bsize);
    printf("Dimensiones de la malla %dx%d\n", Ccolg, Crowg);
    Ngrid = Crowg*Ccolg;
    // Voy a tener como mucho un mínimo por rejilla..
    minGridIDs = malloc(Ngrid*sizeof(int));

    // Para no tener que volver a guardar en memoria y repetir la ejecución
    // TODO: Que vaya cambiando los parametros Wsize, Bsize y Overlap
    clock_t t_init;
    sleep(20);
    while(bucle_entrada){

        printf("\nN cells:    %u\n", Ncells);
        // t_func=omp_get_wtime();
        // t_stage=omp_get_wtime();
        t_init = clock();
        printf("/////////////////////////// LOOP ///////////////////////////\n\n");

        // Me devuelve un mínimo por cada ventana no descartada y guarda el ID en minIDs
        countMin = stage1s(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, octreeIn, min);

        printf("\nCeldas no descartadas:   %d\n", countMin);

        printf("\n\n/////////////////////////// END ///////////////////////////\n");
        // printf("Time elapsed at STAGE 1:     %.6f s\n\n", omp_get_wtime()-t_stage);

        // Para el caso de no hacer solpado; que searcher tenga un valor
        searcher=countMin;

        // Descarto mínimos si hay solape
        // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
        if(Overlap != 0){
            // t_stage=omp_get_wtime();
            printf("\nN minimos de entrada    %d\n", searcher);
            printf("/////////////////////////// MIN SELECT ///////////////////////////\n\n");

            // Ordeno el array de IDs
            qsort(minIDs,countMin,sizeof(int),&cmpfunc);
            // Me quedo solo con los mínimos que se han repetido más de una vez
            searcher = stage2(countMin, minIDs);

            printf("\nNumero de minimos que me quedo: %d \n", searcher);

            printf("\n\n/////////////////////////// END ///////////////////////////\n");
            // printf("Time elapsed at STAGE 2:     %.6f s\n\n",omp_get_wtime() - t_stage );
        }


        // Aplico la malla para ver si hay zonas sin puntos
        if(Bsize > 0){
            printf("N cells grid:    %u\n", Ngrid);
            // t_stage=omp_get_wtime();
            printf("/////////////////////////// GRID ///////////////////////////\n\n\n");

            // Creo un nuevo octree con todos los mínimos; las mismas dimensiones que el grande
            grid = createOctree(center, maxRadius);
            for(int i = 0; i < searcher; i++)
               insertPoint(&pointer[minIDs[i]], grid);

            //
            addMin = stage3s(Bsize, Crowg, Ccolg, minGridIDs, octreeIn, grid, min);

            printf("Minimos añadidos:        %d\n", addMin);
            printf("\n\n/////////////////////////// END ///////////////////////////\n");
            // printf("Time elapsed at STAGE 3:     %.6f s\n\n",omp_get_wtime() - t_stage );
            // Ya no necesito este octree
            free(grid);
            grid = NULL;
        }

        printf("Finalmente, el mapa va a tener %d puntos, %d puntos menos\n", searcher+addMin, Npoints - searcher+addMin );
        printf("TOTAL time elapsed:     %.6f s\n", resultados[--bucle_entrada] = (double)(clock()-t_init) / CLOCKS_PER_SEC );


        if(bucle_entrada){
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
