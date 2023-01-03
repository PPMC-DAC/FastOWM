#define TBB_PREVIEW_GLOBAL_CONTROL 1

#include "../include/envi.h"

#include <chrono>
#include "tbb/global_control.h"


int main( int argc, char* argv[]) {

  //Listas
  Lpoint* point_cloud = NULL;

  unsigned int Npoints=0, Ncells, Ngrid;
  Vector2D center, radius, min, max;
  float maxRadius = 0.0;
  double Width, High, Density, Displace;
  unsigned short  minNumPoints, Crow, Ccol, Crowg, Ccolg;

  unsigned int countMin = 0;

  unsigned int numLLPs = 0;

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
  int NUM_PROCS = (argc>5)? atoi(argv[5]) : 1;
  // Número de repeticiones del experimento
  int numRuns = (argc>6)? atoi(argv[6]) : 1;
  double * results = new double[numRuns];
  // Radio mínimo del nodo hoja
  float minRadius = (argc>7)? atof(argv[7]) : 0.1;

  int MaxNumber = (argc>8)? atoi(argv[8]) : -1;
  
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, NUM_PROCS);

  printf("Input.txt: %s ---> EX. CON %d CORES\n", inputTXT.c_str(), NUM_PROCS);

  if( inputTXT.find("INAER_2011_Alcoy.xyz") != std::string::npos ){ // Alcoy mini
    Npoints = 2772832;
    min.x   = 715244.96;
    max.x   = 716057.75;
    min.y   = 4286623.63;
    max.y   = 4287447.70;
  } else if( inputTXT.find("INAER_2011_Alcoy_Core.xyz") != std::string::npos ){ // Alcoy
    Npoints = 20380212;
    min.x   = 714947.98;
    max.x   = 716361.06;
    min.y   = 4286501.93;
    max.y   = 4288406.23;
  } else if( inputTXT.find("BABCOCK_2017_Arzua_3B.xyz") != std::string::npos ){ //Arzua
    Npoints = 40706503;
    min.x   = 568000.00;
    max.x   = 568999.99;
    min.y   = 4752320.00;
    max.y   = 4753319.99;
  } else if( inputTXT.find("V21_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion forestal
    Npoints = 42384876;
    min.x   = 526964.093;
    max.x   = 527664.647;
    min.y   = 4742610.292;
    max.y   = 4743115.738;
  } else if( inputTXT.find("V19_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion urban
    Npoints = 48024480;
    min.x   = 526955.908;
    max.x   = 527686.445;
    min.y   = 4742586.025;
    max.y   = 4743124.373;
  } else if( inputTXT.find("sample24.xyz") != std::string::npos ){
    Npoints = 7492;
    min.x   = 513748.12;
    max.x   = 513869.97;
    min.y   = 5403124.76;
    max.y   = 5403197.20;
  } else if ( inputTXT.find("Arzua.xyz") == std::string::npos &&
              inputTXT.find("Alcoy.xyz") == std::string::npos && 
              inputTXT.find("BrionF.xyz") == std::string::npos && 
              inputTXT.find("BrionU.xyz") == std::string::npos ){
    printf("No header data!\n");
    exit(-1);
  }

  printf("Reading LiDAR points...\n");
  if(readXYZfile(inputTXT, point_cloud, Npoints, min, max) == -1){
    printf("Unable to read file!\n");
    exit(-1);
  }

  printf("Npoints=%d; xmin = %.2lf; xmax = %.2lf; ymin = %.2lf; ymax = %.2lf\n",Npoints,min.x,max.x,min.y,max.y );

  radius = getRadius(min, max, &maxRadius);
  center = getCenter(min, radius);
  printf("QTREE PARAMETERS:\n");
  printf("MaxRadius:  %.3f\n", maxRadius);
  printf("Center:     %.2f , %.2f\n", center.x,center.y);
  printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
  printf("MinRadius:   %.2f\n", minRadius);
  printf("Building QTREE...\n");

  Width = round2d(max.x-min.x);
  High = round2d(max.y-min.y);
  Density = Npoints/(Width*High);

  // Defines the maximum number of points a leaf node can contain
  if(MaxNumber==-1) MaxNumber = (int)(minRadius*minRadius*Density);

  using tempo_t = std::chrono::steady_clock;
  using cast_t = std::chrono::duration<double, std::milli>;
  std::chrono::time_point<tempo_t> e_tree, s_func, e_stage1, e_stage2, e_stage3;

  printf("INSERTING POINTS...\n");

  e_tree = tempo_t::now();
  Qtree qtreeIn = new Qtree_t( center, maxRadius);
  for(int i = 0; i < Npoints; i++){
    insertPointF2(&point_cloud[i], qtreeIn, minRadius, MaxNumber);
  }
  std::cout << "Time elapsed at Quadtree construction: " << cast_t(tempo_t::now() - e_tree).count() << " ms\n\n";
  
  printf("CLOUD PARAMETERS:\n");
  printf("Number of LiDAR points      %d\n",Npoints);
  printf("Width:  %.2lf\n",Width);
  printf("Height:  %.2lf\n",High);
  printf("Density:  %.3lf\n",Density);

  printf("\nSize of sliding window (SW): %u\n", Wsize);
  printf("Grid size     %u\nOverlap                %.2f\n", Bsize,Overlap);

  minNumPoints = 0.5*Density*Wsize*Wsize;
  printf("Minium number of points in the SW to be considered::   %u\n", minNumPoints);

  Displace = round2d(Wsize*(1-Overlap));
  printf("SW x and y Displacement:  %.2f\n", Displace);

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
  printf("Number of SWs per row:   %d\n",Crow);
  Ncells = Crow*Ccol;
  printf("Total number of OWM steps (Crow x CCol):   %d\n",Ncells);

  // Stage 3 parameters
  printf("\nGrid (for stage 3) parameters:\n");
  Crowg=(int)floor(Width/Bsize)+1;
  Ccolg=(int)floor(High/Bsize)+1;
  printf("Grid dimesions in %dx%d boxes: %dx%d\n\n", Bsize, Bsize, Ccolg, Crowg);
  Ngrid = Crowg*Ccolg;

  // Vectores para los mínimos
  std::vector<int> minIDs;
  std::vector<int> minGridIDs;

  printf("/////////////////////////////////////////////////////////////////////\n");

  while(numRuns){

    s_func = tempo_t::now();

    // stage1s(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min);
    // stage1rem(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min);
    // stage1remCpp(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min);
    // stage1cpp(Wsize, Overlap, Crow, Ccol, minNumPoints, minIDs, qtreeIn, min);
    // minIDs = stage1tbb(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    minIDs = stage1tbbRem(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1reduce(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1reduceRem(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1tg(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1tgRem(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);
    // minIDs = stage1tg2(Wsize, Overlap, Crow, Ccol, minNumPoints, qtreeIn, min);

    e_stage1 = tempo_t::now();

    // Para el caso de no hacer solpado; que numLLPs tenga un valor
    numLLPs = minIDs.size();

    printf("Number of found minima: %d\n", numLLPs);

    if(Overlap != 0.0){

        // Ordeno el array de IDs
        std::sort(minIDs.begin(), minIDs.end(), [](int& a, int& b){
          return a < b;
        });

        // Me quedo solo con los mínimos que se han repetido más de una vez
        numLLPs = stage2(numLLPs, minIDs);

        printf("\nNumber of found LLPs: %d \n", numLLPs);

        e_stage2 = tempo_t::now();
    }

    if(Bsize > 0){

      Qtree grid =  new Qtree_t( center, maxRadius) ;

      if(Overlap != 0.0){
        for(int i = 0; i < numLLPs; i++)
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
      printf("\nNumber of seed points added at stage 3: %d\n\n", addMin);

      // Ya no necesito este qtree
      deleteQtree(grid);
      delete(grid);

    }// Bsize

        

    // printf("REP %d\n", numRuns);
    std::cout << "STAGE 1 time elapsed: " << cast_t(e_stage1 - s_func).count() << " ms\n\n";
    std::cout << "STAGE 2 time elapsed: " << cast_t(e_stage2 - e_stage1).count() << " ms\n\n";
    std::cout << "STAGE 3 time elapsed: " << cast_t(e_stage3 - e_stage2).count() << " ms\n\n";
    std::cout << "TOTAL time elapsed: " << cast_t(e_stage3 - s_func).count() << " ms\n";
    results[--numRuns] = cast_t(e_stage3 - s_func).count()/1e3;

    printf("Output ground seed-point cloud with %d points, %d fewer points than input cloud\n", numLLPs+addMin, Npoints - numLLPs - addMin );

    if(numRuns){
      minIDs.clear();
      minGridIDs.clear();
    }

    printf("/////////////////////////////////////////////////////////////////////\n\n");

  } // while de numRuns

  printf("Time of each run:  ");
  printf("  %.4lf  ", results[0]);
  numRuns = (argc>6)? atoi(argv[6]) : 1;
  double best = results[0];
  if(numRuns > 1){
    for( int i=1 ; i<numRuns ; i++ ){
      printf("  %.4lf  ", results[i]);
      results[0] += results[i];
      if(best > results[i])
        best = results[i];
    }
    printf("\nAverage: %.4lf ms.\t Best time: %4lf ms.\n\n", results[0]/numRuns, best);
  } else 
    printf("\nAverage: %.4lf ms.\t Best time: %4lf ms.\n\n", results[0], best);
    

  // Append vector
  minGridIDs.insert(minGridIDs.end(), minIDs.begin(), minIDs.begin()+numLLPs);
  // check results
  if (check_results(gold_results, minGridIDs, &point_cloud, Displace) < 0) {
      printf("Unable to check results\n");
  }
#ifdef DEBUG
  if(save_file(outputTXT, minGridIDs, &point_cloud) < 0){
    printf("Unable to create file!\n");
  }
#endif

  // Libero memoria
  minIDs.clear();
  minGridIDs.clear();
  free(point_cloud);
  point_cloud = NULL;
  
  deleteQtree(qtreeIn);
  delete(qtreeIn);

  return 0;

}
