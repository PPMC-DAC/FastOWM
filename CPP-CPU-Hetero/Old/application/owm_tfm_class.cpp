#include "../include/envi_class.h"

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int main( int argc, char* argv[]) {

 // Ficheros
    FILE* fileLAS;
    FILE* fileMin;

    // Octrees
    // Octree octreeIn = NULL;
    // Octree grid = NULL;

    // uOctree octreeIn;
    // uOctree grid;

    uOctree grid;

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
    // unsigned int bucle_entrada = 1;

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
    // if(argc>6) bucle_entrada = atoi(argv[6]);
    int bucle_entrada = (argc>6)? atoi(argv[6]) : 1;
    float minRadius = (argc>7)? atof(argv[7]) : 1.0;

    double resultados[bucle_entrada];

    omp_set_num_threads(num_procs);
    // tbb::task_scheduler_init init(num_procs);

    printf("Input.txt: %s ---> EX. CON %d CORES\n", inputTXT, num_procs);

    // Abro el fichero
    if((fileLAS = fopen(inputTXT,"r")) == NULL){
      printf("Unable to open file!\n");
      exit(-1);
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
    } else if(!strcmp(inputTXT,"/mnt/lustre/scratch//home/ulc/cursos/curso357/BABCOCK_2017_Arzua_3B.xyz")){
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
    } else if(!strcmp(inputTXT,"/mnt/lustre/scratch//home/ulc/cursos/curso357/V21_group1_densified_point_cloud.xyz")){
      Npoints = 42384876;
      min.x   = 526964.093;
      max.x   = 527664.647;
      min.y   = 4742610.292;
      max.y   = 4743115.738;
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
    } else {
      printf("No header data\n");
      exit(-1);
    }

    // Reservo memoria para la nube de puntos
    pointer = (Lpoint*)malloc(Npoints*sizeof(Lpoint));

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
        exit(-1);
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
      exit(-1);
    }

    printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );

    // Vector3D mmin, mmax;
    // mmin.x = min.x + minRadius;
    // mmin.y = min.y + minRadius;
    // mmax.x = max.x - minRadius;
    // mmax.y = max.y - minRadius;

    // std::unique_ptr<Octree_t> octetooos[8];
    //
    // octetooos[0] = std::unique_ptr<Octree_t>{(Octree_t*)malloc(sizeof(Octree_t))};
    //
    // octetooos[0]->numPts = 0;

    // nHoja unaHoja;
    // nIntermedio unoInter;

    // std::vector<std::unique_ptr<nIntermedio>> tareas;
    //
    // tareas.push_back(std::unique_ptr<nIntermedio>{ new nHoja });
    // tareas.push_back(std::unique_ptr<nIntermedio>{ new nIntermedio });
    // printf("myfloat en hoja %g\n", tareas[0]->radius);
    // printf("myfloat en hoja %g\n", tareas[1]->radius);
    // tareas[0]->vpoints.push_back(&pointer[0]);
    // uOctree2 octreePrueba( new nIntermedio2(center, maxRadius) );
    // prueba1(octreePrueba);
    // // printf("myfloat %g\n", octreePrueba->radius);
    // // octreePrueba.reset( new nIntermedio(center, maxRadius+1.0) );
    // printf("myfloat %g\n", octreePrueba->radius);
    // std::vector<Lpoint*> v = octreePrueba->getVector();
    // v.push_back(&pointer[0]);
    // printf("myvector size %zu\n", v.size());
    // octreePrueba.reset();

    radius = getRadius(min, max, &maxRadius);
    center = getCenter(min, radius);
    printf("OCTREE PARAMETERS:\n");
    printf("MaxRadius:  %.3f\n", maxRadius);
    printf("Center:     %.2f , %.2f\n", center.x,center.y);
    printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
    printf("InRadius:   %.2f\n", minRadius);
    printf("CREANDO OCTREE...\n");
    // octreeIn = createOctreeF(center, maxRadius);
    // uOctree octreeIn( new nIntermedio(center, maxRadius) );
    uOctree octreeIn( new nHoja2(center, maxRadius) );


    // Inserto los puntos en el Octree
    // for(int i = 0; i < Npoints; i++)
    //    insertPoint(&pointer[i], octreeIn);
    printf("INSERTING POINTS...\n");
    // sleep(1);
    for(int i = 0; i < Npoints; i++){
      // insertPointMinRadius(&pointer[i], octreeIn, minRadius);
      insertPointF(&pointer[i], octreeIn, minRadius, 0);
      // printf("punto %d\n", i);
    }

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
    // minIDs = (int*)malloc(Ncells*sizeof(int));

    // PARAMETROS DE LA MALLA
    printf("\nMALLA:\n");
    Crowg=(int)floor(Width/Bsize)+1;
    Ccolg=(int)floor(High/Bsize)+1;
    // Crow=(int)round(Width/Bsize);
    // Ccol=(int)round(High/Bsize);
    printf("Dimensiones de la malla %dx%d\n", Ccolg, Crowg);
    Ngrid = Crowg*Ccolg;
    // Voy a tener como mucho un mínimo por rejilla..
    // minGridIDs = (int*)malloc(Ngrid*sizeof(int));

    // free(pointer);
    // free(minGridIDs);
    //
    // return 0;

    // Para no tener que volver a guardar en memoria y repetir la ejecución
    // TODO: Que vaya cambiando los parametros Wsize, Bsize y Overlap
    // sleep(5);
    while(bucle_entrada){

        minIDs = (int*)malloc(Ncells*sizeof(int));
        memset(minIDs, -1, Ncells*sizeof(int));
        minGridIDs = (int*)malloc(Ngrid*sizeof(int));
        memset(minGridIDs, -1, Ngrid*sizeof(int));


        // printf("\nN cells:    %u\n", Ncells);
        // Voy a tener como mucho un mínimo por ventana..
        // minIDs = malloc(Ncells*sizeof(int));
        t_func=omp_get_wtime();
        t_stage=omp_get_wtime();
        // printf("/////////////////////////// LOOP ///////////////////////////\n\n");


        // Me devuelve un mínimo por cada ventana no descartada y guarda el ID en minIDs
        // countMin = stage1(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, octreeIn, min);
        printf("STAGE 1 ");
        // sleep(1);
        countMin = stage1cpp(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, octreeIn.get(), min);
        // countMin = stage1tbb(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, octreeIn, min);

        // printf("\nCeldas no descartadas:   %d\n", countMin);

        // printf("\n\n/////////////////////////// END ///////////////////////////\n");
        printf("time elapsed:     %.6f s\n\n", omp_get_wtime()-t_stage);

        // Para el caso de no hacer solpado; que searcher tenga un valor
        searcher=countMin;

        // Descarto mínimos si hay solape
        // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
        if(Overlap != 0.0){
            t_stage=omp_get_wtime();
            // printf("\nN minimos de entrada    %d\n", searcher);
            // printf("/////////////////////////// MIN SELECT ///////////////////////////\n\n");

            // Ordeno el array de IDs
            qsort(minIDs,Ncells,sizeof(int),&cmpfunc);
            // Me quedo solo con los mínimos que se han repetido más de una vez
            // searcher = stage2(countMin, minIDs);
            searcher = stage2(Ncells, minIDs);

            // printf("\nNumero de minimos que me quedo: %d \n", searcher);

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
            // grid = createOctreeF(center, maxRadius);
            uOctree grid( new nHoja2(center, maxRadius) );

            if(Overlap != 0.0){
              for(int i = 0; i < searcher; i++)
                insertPointF(&pointer[minIDs[i]], grid, minRadius, 0);
            } else { // because data is all out of order; NO stage2
              for(int i = 0; i < Ncells; i++)
                if(minIDs[i] >= 0) insertPointF(&pointer[minIDs[i]], grid, minRadius, 0);
            }

            // addMin = stage3s(Bsize, Crowg, Ccolg, minGridIDs, octreeIn, grid, min);
            // addMin = stage3(Bsize, Crowg, Ccolg, minGridIDs, octreeIn, grid, min);
            addMin = stage3cpp(Bsize, Crowg, Ccolg, minGridIDs, octreeIn.get(), grid.get(), min);

            // printf("Minimos añadidos:        %d\n", addMin);
            // printf("\n\n/////////////////////////// END ///////////////////////////\n");
            printf("Time elapsed at STAGE 3:     %.6f s\n\n",omp_get_wtime() - t_stage );
            // Ya no necesito este octree
            deleteOctree(grid.get());
            grid.reset(NULL);
        }

          // printf("REP %d\n", bucle_entrada);
        printf("TOTAL time elapsed:     %.6f s\n", resultados[--bucle_entrada] = omp_get_wtime() - t_func);
        printf("Finalmente, el mapa va a tener %d puntos, %d puntos menos\n", searcher+addMin, Npoints - searcher+addMin );


        if(bucle_entrada){
          free(minIDs);
          minIDs=NULL;
          free(minGridIDs);
          minGridIDs=NULL;
          // Para parar un poco entre vuelta y vuelta
          sleep(5);
        }

        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");
        printf("/////////////////////////////////////////////////////////////////////\n");

    }

    printf("Ejecuciones:  ");
    bucle_entrada = (argc>6)? atoi(argv[6]) : 1;
    for( int i=0 ; i<bucle_entrada ; i++ ){
      printf("  %.4lf  ", resultados[i]);
      if(resultados[0] > resultados[i])
          resultados[0] = resultados[i];
    }
    printf("\nBEST: %.4lf\n", resultados[0]);

    // Fichero de salida
    printf("Creo el fichero %s ...\n", outputTXT);
    if((fileMin = fopen(outputTXT,"w")) == NULL){
      printf("Unable to create file!\n");
      exit(-1);
    }

    if(Overlap != 0.0){
      for(int i=0 ; i<searcher ; i++)
        fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);
    } else { // because data is all out of order; NO stage2
      for(int i=0 ; i<Ncells ; i++)
        if(minIDs[i] >= 0)
          fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minIDs[i]].x, pointer[minIDs[i]].y,pointer[minIDs[i]].z);
    }

    for(int i=0 ; i<addMin ; i++)
      if(minGridIDs[i] >= 0)
        fprintf(fileMin, "%.2f %.2f %.2f\n", pointer[minGridIDs[i]].x, pointer[minGridIDs[i]].y,pointer[minGridIDs[i]].z);

    //Ya no necesito mas el fichero
    if(fclose(fileMin)){
      printf("Cannot close the file\n");
      exit(-1);
    }

    // Libero memoria
    free(minIDs);
    minIDs=NULL;
    free(minGridIDs);
    minGridIDs=NULL;
    free(pointer);
    pointer=NULL;
    deleteOctree(octreeIn.get());
    octreeIn.reset(NULL);

    return 0;
}
