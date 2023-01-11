#include "../src/envi_gpu.hpp"

#ifdef DYNAMIC
#include "../schedulers/LidarSchedulerOneApi.h"
#endif


int main( int argc, const char* argv[]) {

  uint32_t Ncells, Ngrid;
  uint32_t nRows, nColsg, nRowsg;
  Vector2D center, radius, min, max;
  float maxRadius = 0.0;
  double Width, High, Density;

  double best_time = 111111.0;
  double best_rate;

  uint32_t countMin = 0;
  uint32_t numLLPs = 0;
  uint32_t addMin = 0;

  using tempo_t = std::chrono::steady_clock;
  using cast_t = std::chrono::duration<double, std::milli>;
  std::chrono::time_point<tempo_t> s_stage1, e_stage1, s_stage2, e_stage2, s_stage3, e_stage3, s_gpu, e_gpu;
  std::chrono::time_point<tempo_t> i_start, i_end;

  cxxopts::Options options("./parallel [options]", "OWM algorithm to idetify the ground surface from LiDAR data");

  options.add_options()
          ("v,verbose", "Dump input/output data", cxxopts::value<bool>()->default_value("false"))
          ("h,help", "Help message")
          ("i,input", "Input file name without .xyz extension (LiDAR data file)",
            cxxopts::value<std::string>()->default_value("data/INAER_2011_Alcoy"))
          ("W,Wsize", "Sliding Window size", cxxopts::value<uint32_t>()->default_value("10"))
          ("B,Bsize", "Grid size", cxxopts::value<uint32_t>()->default_value("20"))
          ("O,Overlap", "Overlap ratio", cxxopts::value<double>()->default_value("0.80"))
          ("n,npes", "Number of threads", cxxopts::value<int>()->default_value("1"))
          ("l,loop", "Number of runs of the OWM algorithm", cxxopts::value<int>()->default_value("1"))
          ("r,radius", "Value of minRadius ()", cxxopts::value<float>()->default_value("1.0"))
          ("b,balancing", "CPU-GPU partition ratio", cxxopts::value<float>()->default_value("0.5"))
          ("s,size", "Value of maxNumber (max number of points per leaf-node)", cxxopts::value<int>()->default_value("32"))
          ("c,chunk", "GPU chunk", cxxopts::value<uint32_t>()->default_value("1024"))
          ("L,level", "Tree level at which the tree creation becomes parallel", 
            cxxopts::value<int>()->default_value("3"))
          ("d,divide_limit", "Factor limit for the depth of parallelization of the tree creation", 
            cxxopts::value<int>()->default_value("16"));

  auto parameters = options.parse(argc, argv);

  Wsize = parameters["Wsize"].as<uint32_t>();
  uint32_t Bsize = parameters["Bsize"].as<uint32_t>();
  Overlap = parameters["Overlap"].as<double>();
  int npes = parameters["npes"].as<int>();
  int bucle_reps = parameters["loop"].as<int>();
  float minRadius = parameters["radius"].as<float>();
  float rate = parameters["balancing"].as<float>();
  int maxNumber = parameters["size"].as<int>();
  uint32_t chunkGPU = parameters["chunk"].as<uint32_t>();
  int max_level = parameters["level"].as<int>();
  int divide_limit = parameters["divide_limit"].as<int>();

  std::string inputTXT = parameters["input"].as<std::string>() + ".xyz";
  std::string outputTXT = parameters["input"].as<std::string>() + "_salida.xyz";
  std::string gold_results = parameters["input"].as<std::string>() + "_salidaGOLD.xyz";

  std::vector<double> results;
  results.reserve(bucle_reps);

  if (parameters.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
  }


#ifdef DYNAMIC
  LidarSchedulerOneApi lip;
  Params p;
  // p.numcpus = 0; /*de esta forma es SOLO GPU*/
  p.numcpus = 8;
  p.numgpus = 1;
  p.gpuChunk = chunkGPU;
  Dynamic * hs = Dynamic::getInstance(&p);
#endif

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, npes);

  printf("Input.txt: %s ---> EX. with %d CORES\n", inputTXT.c_str(), npes);

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

  } else {
    printf("No header data!\n");
    exit(-1);
  }


  try {
    point_cloud = static_cast<Lpoint*>(mallocWrap(Npoints * sizeof(Lpoint)));
  } catch (sycl::exception &E) {
    std::cout << "point_cloud malloc: " << E.what() << std::endl;
  }

  printf("VOLCANDO PUNTOS...\n");
  // Fichero de entrada
  // if(read_points(inputTXT, &point_cloud) != Npoints){
  if(read_pointsC(inputTXT, point_cloud) < 0){
    printf("Unable to read file!\n");
    // exit(-1);
  }

  printf("xmin = %.2lf\nxmax = %.2lf\nymin = %.2lf\nymax = %.2lf\n",min.x,max.x,min.y,max.y );

  // Dimensiones de la nube para crear el arbol
  radius = getRadius(min, max, &maxRadius);
  center = getCenter(min, radius);
  printf("QTREE PARAMETERS:\n");
  printf("MaxRadius:  %.3f\n", maxRadius);
  printf("Center:     %.2f , %.2f\n", center.x,center.y);
  printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
  printf("minRadius:   %.2f\n", minRadius);
  printf("maxPoints:   %d\n\n", maxNumber);
  // printf("CREANDO QTREE...\n");

  Width = round2d(max.x-min.x);
  High = round2d(max.y-min.y);
  // Densidad en puntos/m^2
  Density = Npoints/(Width*High);

  printf("INSERTING POINTS...\n");
  i_start = tempo_t::now();


  // Qtree cpu_qtree = new Qtree_t( center, maxRadius, NULL );
  // Qtree cpu_qtree = new Qtree_t( NULL, center, maxRadius);
  // cpu_qtree = createQtree(NULL, center, maxRadius);

  // for(int i = 0; i < Npoints; i++){
  //   // insertPoint(&point_cloud[i], cpu_qtree, minRadius);
  //   insertPoint2(&point_cloud[i], cpu_qtree, maxNumber);
  // }

  // cpu_qtree = parallel_qtree_creation(max_level, center, maxRadius, maxNumber);
  cpu_qtree = parallel_qtree_pf2(max_level, center, maxRadius, maxNumber);
  // cpu_qtree = parallel_qtree_pf2(max_level, center, maxRadius, minRadius);
  // cpu_qtree = parallel_qtree_pf3(max_level, center, maxRadius, maxNumber, maxNumber*divide_limit);
  // cpu_qtree = parallel_qtree_creationtg(max_level, center, maxRadius, maxNumber);
  // cpu_qtree = parallel_qtree_creationtg2(center, maxRadius, max_level, maxNumber);

  double cpu_tree_time = cast_t(tempo_t::now() - i_start).count()/1e3;
  
  printf("  CPU Reserved nodes: %zu\n", cpu_tree_nodes.load());

  // i_end = tempo_t::now();
  // std::cout << "INSERT CPU time elapsed: " << cast_t(i_end - i_start).count() << "ms\n";
  std::cout << "  INSERT CPU time elapsed: " << cpu_tree_time << " s\n";

#ifdef INDEX
  try {
    array_indexes = static_cast<QtreeG5>(mallocWrap( cpu_tree_nodes * sizeof(QtreeG5_t)));
  } catch (sycl::exception &E) {
    std::cout << E.what() << std::endl;
  }
#else
  try {
    array_pointers = static_cast<QtreeG4>(mallocWrap( cpu_tree_nodes * sizeof(QtreeG4_t)));
  } catch (sycl::exception &E) {
    std::cout << E.what() << std::endl;
  }
#endif

  try {
    // array_all_points = static_cast<Lpoint**>(mallocWrap( Npoints * sizeof(Lpoint*)));  
    array_all_points = static_cast<Lpoint*>(mallocWrap( Npoints * sizeof(Lpoint)));
  } catch (sycl::exception &E) {
    std::cout << E.what() << std::endl;
  }

  // Qtree parallel_qtree = parallel_qtree_creation(2, center, maxRadius, point_cloud, maxNumber);
  // deleteQtree(parallel_qtree);
  // delete(parallel_qtree);

  // tbb::memory_pool< std::allocator<int> > my_pool;
  // typedef tbb::memory_pool_allocator<int> pool_allocator_t;
  // std::vector<int, pool_allocator_t> myvector( (pool_allocator_t)(my_pool) );

  // tbb::concurrent_vector<int> my_vectors[3];

  // my_vectors[2].push_back(1);

  i_start = tempo_t::now();

#ifdef INDEX
  cudaQtree(cpu_qtree);
#else
  aligned_qtree = launchCopy(cpu_qtree);
#endif

  double new_tree_time = cast_t(tempo_t::now() - i_start).count()/1e3;

  // printf("  Number of points saved in the new Qtree: %d\n", aligned_qtree->numPts);
  // printf("  Number of points saved in the cudaQtree: %d\n", array_indexes[0].numPts);

  std::cout << "  COPY TREE time elapsed: " << new_tree_time << " s\n\n";


  printf("CLOUD PARAMETERS:\n");
  printf("Número de puntos      %lu\n",Npoints);
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
    nCols=(int)(round((Width+2*Wsize*Overlap)/Displace))-1;
    nRows=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
  } else {
    nCols=(int)floor(Width/Wsize)+1;
    nRows=(int)floor(High/Wsize)+1;
  }
  printf("Celdas por columa: %d\n",nRows);
  printf("Celdas por fila:   %d\n",nCols);
  Ncells = nCols*nRows;
  
  // Stage 3 parameters
  printf("\nMALLA:\n");
  nColsg=(int)floor(Width/Bsize)+1;
  nRowsg=(int)floor(High/Bsize)+1;
  printf("Dimensiones de la malla %dx%d\n\n", nRowsg, nColsg);
  Ngrid = nColsg*nRowsg;

  initX = min.x - Wsize/2 + Displace;
  initY = min.y - Wsize/2 + Displace;

  // Vectores para los mínimos
  // std::vector<int> minIDs;
  std::vector<int> minGridIDs;

  printf("/////////////////////////////////////////////////////////////////////\n\n");

  // Para no tener que volver a guardar en memoria y repetir la ejecución
  // TODO: Que vaya cambiando los parametros Wsize, Bsize y Overlap
  // sleep(5);

  // int* minIDs = NULL;
  int countGPU=0;


  while(bucle_reps){

    minIDs = static_cast<int*>(mallocWrap( Ncells*sizeof(int) ));
    memset(minIDs, -1, Ncells*sizeof(int));
    countGPU = 0;

    // s_stage1 = tempo_t::now();

    // Me devuelve un mínimo por cada ventana no descartada y guarda el ID en minIDs

    // stage1s(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min);
    // stage1rem(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, cpu_qtree, min);
    // stage1remCpp(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, cpu_qtree, min);
    // stage1cpp(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, cpu_qtree, min);
    // minIDs = stage1tbb(Wsize, Overlap, nCols, nRows, minNumPoints, cpu_qtree, min);
    // minIDs = stage1tbbRem(Wsize, Overlap, nCols, nRows, minNumPoints, cpu_qtree, min);
    // minIDs = stage1reduce(Wsize, Overlap, nCols, nRows, minNumPoints, cpu_qtree, min);
    // minIDs = stage1reduceRem(Wsize, Overlap, nCols, nRows, minNumPoints, cpu_qtree, min);
    // minIDs = stage1tg(Wsize, Overlap, nCols, nRows, minNumPoints, cpu_qtree, min);
    // minIDs = stage1tgRem(Wsize, Overlap, nCols, nRows, minNumPoints, cpu_qtree, min);
    // minIDs = stage1tg2(Wsize, Overlap, nCols, nRows, minNumPoints, cpu_qtree, min);

    // e_stage1 = tempo_t::now();

    s_gpu = tempo_t::now();

#ifdef DYNAMIC
    hs->heterogeneous_parallel_for(0, Ncells, &lip);
#else
    // stage1s(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min);
    // stage1tbbOne(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min);
    // minIDs = stage1reduce(Wsize, Overlap, nCols, nRows, minNumPoints, aligned_qtree, min);
    // stage1gpu(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, qtreeGPU2, min);
    // stage1gpu(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, qtreeGPU3, min);
    // stage1gpu(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min);
    // stage1gpuRem(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min);
    // stage1heter(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min, rate);
    // stage1hetertg(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min, rate);
    // stage1heterOne(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min, rate);
    // stage1heterTgMix(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, min, rate);
    // stage1heter2q(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min, rate);
    // cudaStage1(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min);
    // stage1index(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, min);
    // stage1heterAcc(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, min, rate);
    // stage1acc(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min, rate);
    stage1tbbRem(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min);
#endif

    e_gpu = tempo_t::now();

    // countMin=0;
    // for(int i=0; i<Ncells; i++){
    //   if(minIDs[i] != -1)
    //     countMin++;
    // }
    // printf("%d celdas; %d minimos no descartados\n", Ncells, countMin);

    // Para el caso de no hacer solpado; que numLLPs tenga un valor
    // numLLPs = minIDs.size();
    // for(int i=0; i<numLLPs; i++){
    //   minGPU[i] = minIDs[i];
    // }

    // Descarto mínimos si hay solape
    // Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
    if(Overlap != 0.0){

        // Me quedo solo con los mínimos que se han repetido más de una vez
        // numLLPs = stage2(numLLPs, minIDs);

        s_stage2 = tempo_t::now();

        qsort(minIDs, Ncells, sizeof(int), &cmpfunc);

        countGPU = stage2GPU(Ncells, minIDs);

        e_stage2 = tempo_t::now();

    }


    // Aplico la malla para ver si hay zonas sin puntos
    if(Bsize > 0){

      s_stage3 = tempo_t::now();

      // Creo un nuevo qtree con todos los mínimos; las mismas dimensiones que el grande
      // Qtree grid =  new Qtree_t( center, maxRadius, NULL) ;
      // Qtree grid =  new Qtree_t( NULL, center, maxRadius) ;
      Qtree grid = createQtree(NULL, center, maxRadius);


      if(Overlap != 0.0){
        for(int i = 0; i < countGPU; i++)
          // insertPoint(&point_cloud[minIDs[i]], grid, 4.0);
          insertPoint2(&point_cloud[minIDs[i]], grid, 128);
      } else { // because data is all out of order; NO stage2
        for(int i = 0; i < Ncells; i++)
          if(minIDs[i] >= 0) insertPoint(&point_cloud[minIDs[i]], grid, 4.0);
      }

      // Qtree grid = parallel_qtree_stage3(max_level, center, maxRadius, 128, minIDs, countGPU);

      // stage3s(Bsize, nColsg, nRowsg, minGridIDs, cpu_qtree, grid, min);
      // stage3cpp(Bsize, nColsg, nRowsg, minGridIDs, cpu_qtree, grid, min);
#ifdef INDEX      
      minGridIDs = stage3reduce(Bsize, nColsg, nRowsg, cpu_qtree, grid, min);
#else
      minGridIDs = stage3reduce(Bsize, nColsg, nRowsg, aligned_qtree, grid, min);
#endif

      e_stage3 = tempo_t::now();

      addMin = minGridIDs.size();
      // printf("\nMinimos añadidos:\t\t%d\n", addMin);

      // Ya no necesito este qtree
      deleteQtree(grid);
      delete(grid);

    }// Bsize

        

    // printf("REP %d\n", bucle_reps);
    // std::cout << "STAGE1 time elapsed: " << cast_t(e_stage1 - s_stage1).count() << "ms\n";
    std::cout << "STAGE1 HETEROGENEOUS " << rate*100 << "% time elapsed: " << cast_t(e_gpu - s_gpu).count() << "ms\n\n";
    std::cout << "STAGE2 time elapsed: " << cast_t(e_stage2 - s_stage2).count() << "ms\n\n";
    std::cout << "STAGE3 time elapsed: " << cast_t(e_stage3 - s_stage3).count() << "ms\n\n";
    double aux_time = cast_t(e_stage3 - s_gpu).count();
    std::cout << "TOTAL time elapsed: " << aux_time << "ms\n";
    // results[--bucle_reps] = aux_time/1e3;
    results.push_back(aux_time/1e3);
    if(aux_time < best_time){
      best_time = aux_time;
      best_rate = rate;
    }

    bucle_reps--;

    // printf("Finalmente, el mapa va a tener %d puntos, %d puntos menos\n", numLLPs+addMin, Npoints - numLLPs+addMin );
    printf("Finalmente, el mapa va a tener %d puntos, %lu puntos menos\n", countGPU+addMin, Npoints - countGPU+addMin );

    if(bucle_reps){

      if(rate > 0.0) rate += 0.01;

      freeWrap(minIDs);
      // freeWrap(minIDs);
      // minIDs = NULL;

      // minIDs.clear();

      minGridIDs.clear();
      // Para parar un poco entre vuelta y vuelta
      sleep(3);
    }

    printf("/////////////////////////////////////////////////////////////////////\n\n");

  } // while de bucle_reps

  printf("Ejecuciones:  ");

  for(double& item : results) printf("  %.4lfs  ", item);
  // printf("\nBEST: %.4lf; minRadius: %g\n", best_time, minRadius);
  printf("\nBEST: %g s (rate = %g); minRadius: %g; maxNumber: %d\n", best_time/1e3, best_rate, minRadius, maxNumber);
  printf("FINAL ( creation + copy + best ): %g s\n", cpu_tree_time + new_tree_time + best_time/1e3);

  // Append vector
  // minGridIDs.insert(minGridIDs.end(), minIDs.begin(), minIDs.begin()+numLLPs);
  for(int i = 0; i < countGPU; i++)
    minGridIDs.push_back(minIDs[i]);
  // check results
  double checked_rate = check_results(gold_results, minGridIDs, &point_cloud, Displace);
  if (checked_rate < 0) {
      printf("Unable to check results\n");
  }

  // Fichero de salida
  if(save_file(outputTXT, minGridIDs, &point_cloud) < 0){
    printf("Unable to create results file!\n");
  }

  if(save_time("results_maxPoints.csv", inputTXT, npes, chunkGPU, minRadius, maxNumber, max_level,
            cpu_tree_time, new_tree_time, best_time/1e3, best_rate, checked_rate) < 0){

    printf("Unable to create results file!\n");

  }

  // Libero memoria
  // minIDs.clear();
  minGridIDs.clear();
  freeWrap(point_cloud);
  // freeWrap(point_cloud);
  // point_cloud = NULL;

  freeWrap(minIDs);
  // freeWrap(minIDs);
  // minIDs = NULL;
  
  deleteQtree(cpu_qtree);
  delete(cpu_qtree);

  // deleteQtreeGPU3(qtreeGPU);
  // delete(qtreeGPU);
  // freeWrap(qtreeGPU,q);

  // printf("Borro el vector global\n");
  // fflush(stdout);
  // all_nodes.clear();
  // all_addresses.clear();

  // for(int i=0; i<cpu_tree_nodes; ++i)
  //   freeWrap(array_pointers[i].quadrants, q);

  // freeWrap(array_all_nodes, q);
#ifdef INDEX
  freeWrap(array_indexes);
#else
  freeWrap(array_pointers);
#endif
  freeWrap(array_all_points);
  // freeWrap(array_pointers);
  // freeWrap(array_indexes);
  // freeWrap(array_all_points);
  aligned_qtree = NULL;

  return 0;

}
