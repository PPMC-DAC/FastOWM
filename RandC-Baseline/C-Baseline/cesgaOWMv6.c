
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

#include "include/environment.h"

double findMinR(Lpoint** A, int* idmin, int n)
{
    double zmin;

    if(n==1)
        return A[0]->z;

    zmin=findMinR(A,idmin,n-1);

    if(zmin<A[n-1]->z)
        return zmin;

    *idmin = n-1;
    return round(A[n-1]->z * 100)/100;
}


int cmpLpoint (const void* a, const void* b){
   // Le pongo el 100 para que se quede con el signo
   // return 100*( round( ((zSearch*)a)->z * 100)/100 - round( ((zSearch*)b)->z * 100)/100 );
   zSearch* punteroA = (zSearch*) a;
   zSearch* punteroB = (zSearch*) b;

   if (round( (punteroB)->z * 100)/100 > round( (punteroA)->z * 100)/100) {
        return -1;
    } else if (round( (punteroB)->z * 100)/100 < round( (punteroA)->z * 100)/100) {
        return 1;
    } else {
        return 0;
    }

}

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

// int compare(const void *_a, const void *_b) {
//
//         int *a, *b;
//
//         a = (int *) _a;
//         b = (int *) _b;
//
//         return (*a - *b);
// }

int main( int argc, char* argv[]){

    // typedef struct{
    //     unsigned int id;
    //     unsigned short count;
    // } accMin;

    // Ficheros
    FILE* fileLAS;
    FILE* fileMin;

    // Octrees
    Octree octreeIn = NULL;
    Octree grid = NULL;

    //Listas
    Lpoint* pointer = NULL;
    int* minIDs = NULL;
    Lpoint** neighbors = NULL;
    // accMin* vMin = NULL;
    Lpoint** minimos = NULL;
    int* minGridIDs = NULL;
    // int* position = NULL;
    //Lista para la salida
    // double* resultados = NULL;

    // zSearch* block_min = NULL;

    unsigned int Npoints;
    Vector3D center, radius, min, max;
    float maxRadius = 0.0;
    double Width, High, Density, Displace;
    unsigned short  minNumPoints, Crow, Ccol;

    Lpoint cellCenter = {0,0.0,0.0,0.0};
    double zzmin;
    unsigned int index=0, countMin=0, cellPoints=0;
    unsigned int idmin = 0;

    double initX, initY;
    int ii,jj; //para el control de los bucles paralelizados

    unsigned int searcher=0;

    unsigned int addMin=0;

    double t_stage, t_func;

    // Tamaño de la ventana deslizante
    unsigned short Wsize = 12;
    // Tamaño de la rejilla
    unsigned short Bsize = 20;
    // Solape de la ventana deslizante
    double Overlap = 0.8;
    // Numero de procesadores
    unsigned short num_procs = 4;
    // Tamaño del bloque del scheduler de OMP
    // unsigned short block_size = 1;

    //Control del bucle de ejecución y parámetros
    unsigned int bucle_entrada = 1;

    char inputTXT[50] = {"../datos/INAER_2011_Alcoy.xyz"};
    // char inputLAS[50] = {"../datos/INAER_2011_Alcoy.las"};
    char outputTXT[50] = {"../datos/INAER_2011_Alcoy_salida.xyz"};
    // char outputLAS[50] = {"../datos/INAER_2011_Alcoy_salida.las"};

    // Compruebo los argumentos de entrada
    if(argc>1) {
      strcpy(inputTXT,argv[1]);
      // strcpy(inputLAS,argv[1]);
      strcat(inputTXT,".xyz");
      // strcat(inputLAS,".las");

      strcpy(outputTXT,argv[1]);
      // strcpy(outputLAS,argv[1]);
      strcat(outputTXT,"_salida.xyz");
      // strcat(outputLAS,"_salida.las");
    }
    if(argc>2) Wsize = atoi(argv[2]);
    if(argc>3) Bsize = atoi(argv[3]);
    if(argc>4) Overlap = atof(argv[4]);
    if(argc>5) num_procs = atoi(argv[5]);
    if(argc>6) bucle_entrada = atoi(argv[6]);

    double resultados[bucle_entrada];

    omp_set_num_threads(num_procs);

    printf("Input.txt: %s ---> EX. CON %d CORES\n", inputTXT, num_procs);
    // printf("Input.las: %s\n", inputLAS);

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


    // Abro el fichero
    if((fileLAS = fopen(inputTXT,"r")) == NULL){
    // if((fileLAS = fopen("../datos/INAER_2011_Alcoy.xyz","r")) == NULL){
    // if((fileLAS = fopen("../datos/INAER_2011_Alcoy_Core.xyz","r")) == NULL){
    // if((fileLAS = fopen("../datos/INSITU_2018_Jenaro.xyz","r")) == NULL){
      printf("Unable to open file!\n");
      return -1;
    }

    // Reservo memoria para la nube de puntos
    pointer = malloc(Npoints*sizeof(Lpoint));

    printf("VOLCANDO PUNTOS...\n");
    for(int i=0; i<Npoints ; i++){
      //Obtengo los datos id X Y Z
      pointer[i].id = i;
      if(fscanf(fileLAS, "%lf%lf%lf",&pointer[i].x,&pointer[i].y,&pointer[i].z) < 3){
        printf("Imposible to obtain values\n");
        return -1;
      }
      //Los siguientes 6 campos no los necesito
      while(fgetc(fileLAS)!='\n');
     	// printf("[%u]  %.2lf, %.2lf, %.2lf\n",pointer[i].id,pointer[i].x,pointer[i].y,pointer[i].z);
    }

    // qsort(pointer,Npoints,sizeof(Lpoint),&cmpLpoint);
    // for(int i=0; i<Npoints; i++) printf(" %.2f ", pointer[i].z); printf("\n\n");

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

    // #pragma omp parallel for firstprivate(Npoints,pointer,octreeIn) schedule(static,1)
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

    // printf("BUCLE? ");
    // scanf("%d", &bucle_entrada);

    while(bucle_entrada){

        // printf("BUCLE?\n");
        // scanf("%d", &bucle_entrada);

        printf("\nTamaño de ventana     %u\n", Wsize);
        printf("Tamaño de rejilla     %u\nSolape                %.2f\n", Bsize,Overlap);

        // El numero minimo sera la mitad del numero de puntos medio por celda
        minNumPoints = 0.5*Density*Wsize*Wsize;
        printf("Minimo numero de puntos por celda:   %u\n", minNumPoints);

        Displace = round2d(Wsize*(1-Overlap));
        printf("Displacement   %.2f\n", Displace);
        // unsigned short NumDisplace = Wsize/Displace-1;
        // printf("N displacements   %u\n", NumDisplace);

        // Crow=(int)round(Width/Wsize);
        // Ccol=(int)round(High/Wsize);
        Crow=(int)floor(Width/Wsize)+1;
        Ccol=(int)floor(High/Wsize)+1;
        if(Overlap!=0) Crow=(int)(round((Width+2*Wsize*Overlap)/Displace))-1, Ccol=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
        printf("Celdas por columa: %d\n",Ccol);
        printf("Celdas por fila:   %d\n",Crow);

        printf("\nN cells:    %u\n", Crow*Ccol);

        // Como voy a tener como mucho un mínimo por celda...
        minIDs = malloc(Crow*Ccol*sizeof(int));

        t_func=omp_get_wtime();
        t_stage=omp_get_wtime();
        printf("/////////////////////////// LOOP ///////////////////////////\n\n");

        initX = min.x + Wsize/2 - Wsize*Overlap;
        initY = min.y + Wsize/2 - Wsize*Overlap;

        // printf("%d\n", omp_get_num_threads());
        // countMin = stage1(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, octreeIn, min);

        countMin=0;
        // #pragma omp parallel for shared(countMin) firstprivate(minIDs,octreeIn,Wsize,Overlap,Displace,initX,initY) \
        //                           private(cellCenter,neighbors,cellPoints,idmin) schedule(dynamic,block_size)
        #pragma omp parallel shared(countMin) firstprivate(Ccol,Crow,minIDs,octreeIn,Wsize,Displace,initX,initY) \
                                  private(ii,jj,cellCenter,neighbors,cellPoints,idmin)
        {
            // printf("Thread: %d\n", omp_get_thread_num());
            for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=num_procs ){

                cellCenter.y = initY + jj*Displace;

                // #pragma omp parallel shared(countMin,minIDs) firstprivate(Crow,cellCenter,octreeIn,Wsize,Overlap,Displace,initX) \
                //                           private(ii,neighbors,cellPoints,idmin) num_threads(1)
                // {
                    // for( ii = omp_get_thread_num() ; ii < Crow ; ii+=1 ){
                    #pragma omp parallel for shared(countMin) schedule(dynamic,1)
                    for( ii=0 ; ii < Crow ; ii++ ){

                        cellCenter.x = initX + ii*Displace;
                        // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
                        // printf("Busco los vecinos\n");
                        neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
                        // neighbors = searchNeighbors2Dmin(&cellCenter, octreeIn, Wsize/2, &cellPoints, &idmin);
                        // printf("Numero de elementos de la celda: %d\n", cellPoints );
                        if(cellPoints >= minNumPoints ){

                            // printf("Numero de elementos de la celda: %d\n", cellPoints );
                            idmin = findMin(neighbors, cellPoints);
                            // idmin=0;
                            // findMinR(neighbors,&idmin,cellPoints);

                            // printf("ID min: %d\n", idmin);

                            #pragma omp critical
                            {
                                minIDs[countMin] = neighbors[idmin]->id;
                                countMin++;
                            }
                            // printf("El mínimo %d de la celda es: %.2f\n",countMin ,neighbors[idmin]->z);
                            // free(block_min);
                            // block_min = NULL;
                        }

                        free(neighbors);
                        neighbors = NULL;
                    }
                // }
            }
        }

        printf("\nCeldas no descartadas:   %d\n", countMin);
        // printf("Minimos seleccionados:\n");
        // for(int i=0; i<countMin; i++) printf(" %.2f ", pointer[minIDs[i]].z); printf("\n");
        // for(int i=0; i<countMin; i++) printf(" %d ", minIDs[i]); printf("\n");

        printf("\n\n/////////////////////////// END ///////////////////////////\n");
        printf("Time elapsed at STAGE 1:     %.6f s\n\n", omp_get_wtime()-t_stage);

        // Para el caso de no hacer solpado, que searcher tenga un valor
        searcher=countMin;

        // Descarto mínimos si hay solape
        // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
        if(Overlap != 0){
            t_stage=omp_get_wtime();
            printf("\nN minimos de entrada    %d\n", searcher);
            printf("/////////////////////////// MIN SELECT ///////////////////////////\n\n");


            qsort(minIDs,countMin,sizeof(int),&cmpfunc);

            index=0;
            for( ii=0 ; ii<countMin ; ii=jj ){

                for( jj=ii+1 ; minIDs[ii]==minIDs[jj] && jj<countMin ; jj++ );

                if(jj-ii > 1){
                  minIDs[index]=minIDs[ii];
                  index++;
                }

            }
            searcher=index;


            /*NO VOY A SALIRME DEL ARRAY PORQUE OCUPA MÁS DE countMin*/
            // printf("modulo: %d\n", 6%5);
            // int rango_trabajo = countMin/num_procs;
            // // Para quedarme siempre dentro
            // if(countMin%num_procs) rango_trabajo = floor(countMin/num_procs);

            // index=0;
            // #pragma omp parallel shared(index) firstprivate(rango_trabajo) private(ii,jj)
            // {
            //   int private_index=0;
            //   int private_min[rango_trabajo];
            //   int inicio = rango_trabajo*omp_get_thread_num();
            //   int final = inicio + rango_trabajo;
            //   for( ii=inicio ; ii<final ; ii=jj ){
            //       // position[index] = ii;
            //       // vMin[index].id = minIDs[ii];
            //       // vMin[index].count = 1;
            //       // printf(" /////NUEVO  :%d ", minIDs[ii]);
            //       for( jj=ii+1 ; minIDs[ii]==minIDs[jj] && jj<final ; jj++ ){
            //           // vMin[index].count++;
            //           // printf(" jj:%d", minIDs[jj]);
            //       }
            //       if(jj-ii-1 > 0){
            //         private_min[private_index]=minIDs[ii];
            //         private_index++;
            //       }
            //       // index++;
            //   }
            //   for( ii=0 ; ii<private_index ; ii++){
            //     #pragma omp critical
            //     {
            //       minIDs[index] = private_min[ii];
            //       index++;
            //     }
            //   }
            // }
            // searcher=index;


            // printf("Finalmente me quedo con: \n");
            // for(int i=0; i<countMin; i++) printf(" %.2f ", pointer[minIDs[i]].z); printf("\n");
            // for(int i=0; i<searcher ; i++) printf(" %d ", minIDs[i]); printf("\n");

            printf("\nNumero de minimos que me quedo: %d ; %d ; %d\n", searcher, index, countMin);

            // free(vMin);
            // vMin = NULL;
            // free(position);
            // position = NULL;
            printf("\n\n/////////////////////////// END ///////////////////////////\n");
            printf("Time elapsed at STAGE 2:     %.6f s\n\n",omp_get_wtime() - t_stage );
        }


        //Aquí es donde tengo que tomar la decisión de si voy con todo o desprecio el ultimo cacho de nube
        Crow=(int)floor(Width/Bsize)+1;
        Ccol=(int)floor(High/Bsize)+1;
        // Crow=(int)round(Width/Bsize);
        // Ccol=(int)round(High/Bsize);
        printf("Dimensiones de la malla %dx%d\n", Ccol,Crow);
        printf("N cells grid:    %u\n", Crow*Ccol);
        // int minGridIDs[Crow*Ccol];
        minGridIDs = malloc(Crow*Ccol*sizeof(int));
        // Aplico la malla para ver si hay zonas sin puntos
        if(Bsize > 0){
            t_stage=omp_get_wtime();
            printf("/////////////////////////// GRID ///////////////////////////\n\n\n");

            //Creo un nuevo octree con todos los mínimos; las mismas dimensiones que el grande
            radius = getRadius(min, max, &maxRadius);
            center = getCenter(min, radius);
            // printf("Creo el octree\n");
            grid = createOctree(center, maxRadius);
            for(int i = 0; i < searcher; i++)
               insertPoint(&pointer[minIDs[i]], grid);
            // printf("Termino de insertar los puntos\n");

            // addMin = stage3(Bsize, Crow, Ccol, minGridIDs, octreeIn, grid, min);

            addMin=0;
            #pragma omp parallel shared(addMin) firstprivate(Ccol,Crow,minGridIDs,octreeIn,grid,Bsize) \
                                      private(ii,jj,cellCenter,neighbors,minimos,cellPoints,idmin)
            {

                for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=num_procs ){
                // for( jj = 0 ; jj < Ccol ; jj++ ){
                   cellCenter.y = min.y + Bsize/2 + jj*Bsize;

                   for( ii = 0 ; ii < Crow ; ii++ ){
                       cellCenter.x = min.x + Bsize/2 + ii*Bsize;

                       minimos = searchNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
                       // printf("Numero de elementos de la celda: %d\n", cellPoints );
                       //Tengo que hacerlo porque el algoritmo no lo hace
                       if(cellPoints == 0){
                           // printf("    Tengo una celda vacia en la malla\n" );
                           neighbors = searchNeighbors2D(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                           if(cellPoints>0){

                               idmin = findMin(neighbors, cellPoints);
                               #pragma omp critical
                               {
                                  minGridIDs[addMin] = neighbors[idmin]->id;
                                  addMin++;
                               }
                               // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
                           }
                           free(neighbors);
                           neighbors = NULL;

                       }
                       free(minimos);
                       minimos=NULL;
                   }
                }
            }
            printf("Minimos añadidos:        %d\n", addMin);
            // for(int i=0; i<addMin ; i++)
            //   printf("%d: %.2f , %.2f , %.2f\n", i+1, pointer[minGridIDs[i]].x, pointer[minGridIDs[i]].y,pointer[minGridIDs[i]].z);
            printf("\n\n/////////////////////////// END ///////////////////////////\n");
            printf("Time elapsed at STAGE 3:     %.6f s\n\n",omp_get_wtime() - t_stage );
            // Ya no necesito este octree
            free(grid);
            grid = NULL;
        }
        // printf("Minimos añad.    %d\n", addMin);

        printf("Finalmente, el mapa va a tener %d puntos, %d puntos menos\n", searcher+addMin, Npoints - searcher+addMin );
        printf("TOTAL time elapsed:     %.6f s\n", resultados[--bucle_entrada] = omp_get_wtime() - t_func);
        // Ya no necesito el octreeIn
        // free(octreeIn);
        // octreeIn=NULL;

        // Si voy a seguir, los libero; de no seguir, los dejo para hacer el fichero
        if(bucle_entrada){
          free(minIDs);
          minIDs=NULL;
          free(minGridIDs);
          minGridIDs=NULL;
          // Espero para la siguiente vuelta
          sleep(10);
        }
        // printf("signo: %d\n", signbit(-0.1*100));
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
      fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minIDs[i]].x, pointer[minIDs[i]].y, pointer[minIDs[i]].z);
    // printf("Minimos etapa3...\n");
    for(int i=0 ; i<addMin ; i++)
      fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minGridIDs[i]].x, pointer[minGridIDs[i]].y, pointer[minGridIDs[i]].z);

    //Ya no necesito mas el fichero
    if(fclose(fileMin)){
      printf("Cannot close the file\n");
      return -1;
    }

    // printf("FIN!\n");

    // Libero memoria
    free(octreeIn);
    octreeIn=NULL;
    free(minIDs);
    minIDs=NULL;
    free(minGridIDs);
    minGridIDs=NULL;
    free(pointer);
    pointer=NULL;
    // free(resultados);
    // resultados=NULL;


    // char command[200];
    // printf("Creo el fichero %s ...\n", outputLAS);
    // printf("txt2las -i %s -o %s --parse xyz\n",outputTXT,outputLAS);
    // printf("> source(\"/home/felipe/Escritorio/Beca_CiTIUS/Proyecto_LiDAR/algoritmo_OWM_LiDAR/pintaLAS.R.R\")\n");
    // sprintf(command,"txt2las -i %s -o %s --parse xyz",outputTXT,outputLAS);
    // // sprintf(command,"python3 script_surface.py %s %s",outputTXT,outputLAS);
    // if(system(command) < 0){
    //   perror("system");
    //   return 1;
    // }


    return 0;
}
