#define TBB_PREVIEW_GLOBAL_CONTROL 1

#include "../include/envi.h"

#include <chrono>
#include "tbb/global_control.h"

double Displace;
unsigned short minNumPoints, Crow;
unsigned short Wsize;
Qtree qtreeIn = NULL;

int main( int argc, char* argv[]) {

  // Ficheros
  FILE* fileLAS;
  FILE* fileMin;

  //Listas
  Lpoint* point_cloud = NULL;

  unsigned int Npoints=0, Ncells, Ngrid;
  Vector2D center, radius, min, max;
  float maxRadius = 0.0;
  double Width, High, Density, Displace;
  unsigned short  minNumPoints, Crow, Ccol, Crowg, Ccolg;

  unsigned int countMin = 0;

  unsigned int searcher = 0;

  unsigned int addMin = 0;

  std::string inputTXT = {"./data/INAER_2011_Alcoy.xyz"};
  std::string outputTXT = {"./data/INAER_2011_Alcoy_salida.xyz"};
  std::string gold_results = {"./data/INAER_2011_Alcoy_salidaGOLD.xyz"};

  // Compruebo los argumentos de entrada
  if(argc>1) {
    inputTXT = argv[1];
    inputTXT += ".xyz";

    outputTXT = argv[1];
    outputTXT += "_salida.xyz";

    gold_results = argv[1];
    gold_results += "_salidaGOLD.xyz";
  }
  // Tamaño de la ventana deslizante
  unsigned short Wsize = (argc>2)? atoi(argv[2]) : 10;
  // Tamaño de la rejilla
  unsigned short Bsize = (argc>3)? atoi(argv[3]) : 20;
  // Solape de la ventana deslizante
  double Overlap = (argc>4)? atof(argv[4]) : 0.8;
  // Número de cores
  int NUM_PROCS = (argc>5)? atoi(argv[5]) : 4;
  // Número de repeticiones del experimento
  int bucle_entrada = (argc>6)? atoi(argv[6]) : 1;
  double resultados[bucle_entrada];
  // Radio mínimo del nodo hoja
  float minRadius = (argc>7)? atof(argv[7]) : 1.0;

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, NUM_PROCS);

  // omp_set_num_threads(NUM_PROCS);

  printf("Input.txt: %s ---> EX. CON %d CORES\n", inputTXT.c_str(), NUM_PROCS);


  if( inputTXT.find("INAER_2011_Alcoy.xyz") != std::string::npos ){ // Alcoy mini
    Npoints = 2772832;
    min.x   = 715244.96;
    max.x   = 716057.75;
    min.y   = 4286623.63;
    max.y   = 4287447.70;
    // min.z   = 0;
    // max.z   = 0; //No lo consulto nunca
  } else if( inputTXT.find("INAER_2011_Alcoy_Core.xyz") != std::string::npos ){ // Alcoy
    Npoints = 20380212;
    min.x   = 714947.98;
    max.x   = 716361.06;
    min.y   = 4286501.93;
    max.y   = 4288406.23;
    // min.z   = 0;
    // max.z   = 0; //No lo consulto nunca
  } else if( inputTXT.find("BABCOCK_2017_Arzua_3B.xyz") != std::string::npos ){ //Arzua
    Npoints = 40706503;
    min.x   = 568000.00;
    max.x   = 568999.99;
    min.y   = 4752320.00;
    max.y   = 4753319.99;
    // min.z   = 0;
    // max.z   = 0; //No lo consulto nunca
  } else if( inputTXT.find("V21_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion forestal
    Npoints = 42384876;
    min.x   = 526964.093;
    max.x   = 527664.647;
    min.y   = 4742610.292;
    max.y   = 4743115.738;
    // min.z   = 0;
    // max.z   = 0; //No lo consulto nunca
  } else if( inputTXT.find("V19_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion urban
    Npoints = 48024480;
    min.x   = 526955.908;
    max.x   = 527686.445;
    min.y   = 4742586.025;
    max.y   = 4743124.373;
    // min.z   = 0;
    // max.z   = 0; //No lo consulto nunca
  } else if( inputTXT.find("sample24.xyz") != std::string::npos ){
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
  point_cloud = (Lpoint*)malloc(Npoints*sizeof(Lpoint));

  printf("VOLCANDO PUNTOS...\n");
  // Fichero de entrada
  if(read_points(inputTXT, &point_cloud) != Npoints){
    printf("Unable to read file!\n");
  }

  printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );

  // Dimensiones de la nube para crear el arbol
  radius = getRadius(min, max, &maxRadius);
  center = getCenter(min, radius);
  printf("QTREE PARAMETERS:\n");
  printf("MaxRadius:  %.3f\n", maxRadius);
  printf("Center:     %.2f , %.2f\n", center.x,center.y);
  printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
  printf("InRadius:   %.2f\n", minRadius);
  printf("CREANDO QTREE...\n");

  Width = round2d(max.x-min.x);
  High = round2d(max.y-min.y);
  // Densidad en puntos/m^2
  Density = Npoints/(Width*High);

  // Defines a maximum number of points from which to divide the cell
  int medSize = (int)(minRadius*minRadius*Density);

  Qtree qtreeIn = new Qtree_t( center, maxRadius);
  // Qtree qtreeIn = new Qtree_t( NULL, center, maxRadius);

  printf("INSERTING POINTS...\n");
  for(int i = 0; i < Npoints; i++){
    // insertPointF(&point_cloud[i], qtreeIn, minRadius);
    insertPointF2(&point_cloud[i], qtreeIn, minRadius, medSize);
  }

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
  
  // Stage 3 parameters
  printf("\nMALLA:\n");
  Crowg=(int)floor(Width/Bsize)+1;
  Ccolg=(int)floor(High/Bsize)+1;
  printf("Dimensiones de la malla %dx%d\n\n", Ccolg, Crowg);
  Ngrid = Crowg*Ccolg;

  // Vectores para los mínimos
  std::vector<int> minIDs;
  std::vector<int> minGridIDs;

  printf("/////////////////////////////////////////////////////////////////////\n");
  printf("/////////////////////////////////////////////////////////////////////\n");
  printf("/////////////////////////////////////////////////////////////////////\n\n");

  // Para no tener que volver a guardar en memoria y repetir la ejecución
  // TODO: Que vaya cambiando los parametros Wsize, Bsize y Overlap
  // sleep(5);
  using tempo_t = std::chrono::steady_clock;
  using cast_t = std::chrono::duration<double, std::milli>;
  std::chrono::time_point<tempo_t> s_func, e_stage1, e_stage2, e_stage3;


  while(bucle_entrada){

    s_func = tempo_t::now();

    // Me devuelve un mínimo por cada ventana no descartada y guarda el ID en minIDs

    // stage1s(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min);
    // stage1rem(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min);
    // stage1remCpp(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min);
    // stage1cpp(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min);
    // minIDs = stage1tbb(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1tbbRem(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    minIDs = stage1reduce(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1reduceRem(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1tg(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1tgRem(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1tg2(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);

    e_stage1 = tempo_t::now();

    // Para el caso de no hacer solpado; que searcher tenga un valor
    searcher = minIDs.size();

    printf("\nCeldas no descartadas:\t\t%d\n", searcher);

    // Descarto mínimos si hay solape
    // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
    if(Overlap != 0.0){

        // Ordeno el array de IDs
        std::sort(minIDs.begin(), minIDs.end(), [](int& a, int& b){
          return a < b;
        });

        // Me quedo solo con los mínimos que se han repetido más de una vez
        searcher = stage2(searcher, minIDs);

        printf("\nNumero de minimos que me quedo:\t%d \n", searcher);

        e_stage2 = tempo_t::now();
    }


    // Aplico la malla para ver si hay zonas sin puntos
    if(Bsize > 0){

      // Creo un nuevo qtree con todos los mínimos; las mismas dimensiones que el grande
      Qtree grid =  new Qtree_t( center, maxRadius) ;
      // Qtree grid =  new Qtree_t( NULL, center, maxRadius) ;

      if(Overlap != 0.0){
        for(int i = 0; i < searcher; i++)
          insertPointF(&point_cloud[minIDs[i]], grid, 4.0);
      } else { // because data is all out of order; NO stage2
        for(int i = 0; i < Ncells; i++)
          if(minIDs[i] >= 0) insertPointF(&point_cloud[minIDs[i]], grid, minRadius);
      }

      // stage3s(Bsize, Crowg, Ccolg, minGridIDs, qtreeIn, grid, min);
      // stage3cpp(Bsize, Crowg, Ccolg, minGridIDs, qtreeIn, grid, min);
      minGridIDs = stage3reduce(Bsize, Crowg, Ccolg, qtreeIn, grid, min);

      e_stage3 = tempo_t::now();

      addMin = minGridIDs.size();
      printf("\nMinimos añadidos:\t\t%d\n", addMin);

      // Ya no necesito este qtree
      deleteQtree(grid);
      delete(grid);

    }// Bsize

        

    // printf("REP %d\n", bucle_entrada);
    std::cout << "STAGE1 time elapsed: " << cast_t(e_stage1 - s_func).count() << "ms\n\n";
    std::cout << "STAGE2 time elapsed: " << cast_t(e_stage2 - e_stage1).count() << "ms\n\n";
    std::cout << "STAGE3 time elapsed: " << cast_t(e_stage3 - e_stage2).count() << "ms\n\n";
    std::cout << "TOTAL time elapsed: " << cast_t(e_stage3 - s_func).count() << "ms\n";
    resultados[--bucle_entrada] = cast_t(e_stage3 - s_func).count()/1e3;

    printf("Finalmente, el mapa va a tener %d puntos, %d puntos menos\n", searcher+addMin, Npoints - searcher+addMin );

    if(bucle_entrada){

      minIDs.clear();

      minGridIDs.clear();
      // Para parar un poco entre vuelta y vuelta
      sleep(5);
    }

    printf("/////////////////////////////////////////////////////////////////////\n");
    printf("/////////////////////////////////////////////////////////////////////\n");
    printf("/////////////////////////////////////////////////////////////////////\n\n");

  } // while de bucle_entrada

  printf("Ejecuciones:  ");
  bucle_entrada = (argc>6)? atoi(argv[6]) : 1;
  for( int i=0 ; i<bucle_entrada ; i++ ){
    printf("  %.4lfs  ", resultados[i]);
    if(resultados[0] > resultados[i])
        resultados[0] = resultados[i];
  }
  printf("\nBEST: %.4lf\n", resultados[0]);

  // Append vector
  minGridIDs.insert(minGridIDs.end(), minIDs.begin(), minIDs.begin()+searcher);
  // check results
  if (check_results(gold_results, minGridIDs, &point_cloud, Displace) < 0) {
      printf("Unable to check results\n");
  }

  // Fichero de salida
  if(save_file(outputTXT, minGridIDs, &point_cloud) < 0){
    printf("Unable to create file!\n");
  }

  // Libero memoria
  minIDs.clear();
  minGridIDs.clear();
  free(point_cloud);
  point_cloud = NULL;
  
  deleteQtree(qtreeIn);
  delete(qtreeIn);

  return 0;

}
