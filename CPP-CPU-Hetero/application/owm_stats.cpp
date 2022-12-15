#include "../src/envi_gpu.hpp"

#include "../include/cxxopts.hpp"
#include <chrono>

int main( int argc, const char* argv[]) 
{

  //Listas
  // Lpoint* point_cloud = NULL;

  // unsigned int Npoints=0, Ncells, Ngrid;
  unsigned int Ncells, Ngrid;
  Vector2D center, radius, min, max;
  float maxRadius = 0.0;
  double Width, High, Density, Displace;
  unsigned short  minNumPoints, Crow, Ccol, Crowg, Ccolg;

  unsigned int countMin = 0;

  unsigned int searcher = 0;

  unsigned int addMin = 0;

  using tempo_t = std::chrono::steady_clock;
  using cast_t = std::chrono::duration<double, std::milli>;
  std::chrono::time_point<tempo_t> s_stage1, e_stage1, s_stage2, e_stage2, s_stage3, e_stage3, s_gpu, e_gpu;
  std::chrono::time_point<tempo_t> i_start, i_end;

  // std::string inputTXT = {"./data/INAER_2011_Alcoy.xyz"};
  // std::string outputTXT = {"./data/INAER_2011_Alcoy_salida.xyz"};
  // std::string gold_results = {"./data/INAER_2011_Alcoy_salidaGOLD.xyz"};

  // // Compruebo los argumentos de entrada
  // if(argc>1) {
  //   inputTXT = argv[1];
  //   inputTXT += ".xyz";

  //   outputTXT = argv[1];
  //   outputTXT += "_salida.xyz";

  //   gold_results = argv[1];
  //   gold_results += "_salidaGOLD.xyz";
  // }
  // // Tamaño de la ventana deslizante
  // unsigned short Wsize = (argc>2)? atoi(argv[2]) : 12;
  // // Tamaño de la rejilla
  // unsigned short Bsize = (argc>3)? atoi(argv[3]) : 20;
  // // Solape de la ventana deslizante
  // double Overlap = (argc>4)? atof(argv[4]) : 0.5;
  // // Número de cores
  // int NUM_PROCS = (argc>5)? atoi(argv[5]) : 2;
  // // Número de repeticiones del experimento
  // int bucle_entrada = (argc>6)? atoi(argv[6]) : 1;
  // double resultados[bucle_entrada];
  // // Radio mínimo del nodo hoja
  // float minRadius = (argc>7)? atof(argv[7]) : 1.0;
  // // Factor de reparto de trabajo GPU-CPU
  // float rate = (argc>8)? atof(argv[8]) : 0.5;
  // // Defines a maximum number of points from which to divide the cell
  // int maxSize = (argc>9)? atoi(argv[9]) : 16;

  cxxopts::Options options("./reductionts [options]", "Application to reduce points based on Visvalingam algorithm");

  options.add_options()
          ("v,verbose", "Show data input/output", cxxopts::value<bool>()->default_value("false"))
          ("h,help", "Displays code help")
          ("i,input", "Input file name without extension where points are saved",
            cxxopts::value<std::string>()->default_value("data/INAER_2011_Alcoy"))
          ("W,Wsize", "Window size", cxxopts::value<unsigned short>()->default_value("12"))
          ("B,Bsize", "Grid size", cxxopts::value<unsigned short>()->default_value("20"))
          ("O,Overlap", "Overlap rate", cxxopts::value<double>()->default_value("0.5"))
          ("n,npes", "Number of cores", cxxopts::value<int>()->default_value("2"))
          ("l,loop", "Number of repetitions", cxxopts::value<int>()->default_value("2"))
          ("r,radius", "Leaf nodes radius value", cxxopts::value<float>()->default_value("1.0"))
          ("b,balancing", "Balancing rate", cxxopts::value<float>()->default_value("0.5"))
          ("s,size", "Nodes maximum point size", cxxopts::value<int>()->default_value("32"))
          ("L,level", "Maximum level up to which to distribute the creation", cxxopts::value<int>()->default_value("3"));

  auto parameters = options.parse(argc, argv);

  std::string inputTXT = parameters["input"].as<std::string>() + ".xyz";
  std::string outputTXT = parameters["input"].as<std::string>() + "_salida.xyz";
  std::string gold_results = parameters["input"].as<std::string>() + "_salidaGOLD.xyz";

  unsigned short Wsize = parameters["Wsize"].as<unsigned short>();
  unsigned short Bsize = parameters["Bsize"].as<unsigned short>();
  double Overlap = parameters["Overlap"].as<double>();
  int NUM_PROCS = parameters["npes"].as<int>();
  int bucle_entrada = parameters["loop"].as<int>();
  float minRadius = parameters["radius"].as<float>();
  float rate = parameters["balancing"].as<float>();
  int maxSize = parameters["size"].as<int>();
  int max_level = parameters["level"].as<int>();

  double resultados[bucle_entrada];

  if (parameters.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
  }

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


  try {
    point_cloud = static_cast<Lpoint*>(mallocWrap(Npoints * sizeof(Lpoint)));
  } catch (cl::sycl::invalid_parameter_error &E) {
    std::cout << E.what() << std::endl;
  }

  // std::atomic<int> at_cpu_tree_nodes = 0;


  printf("VOLCANDO PUNTOS...\n");
  // Fichero de entrada
  // if(read_points(inputTXT, &point_cloud) != Npoints){
  //   printf("Unable to read file!\n");
  // }
  if(read_pointsC(inputTXT, point_cloud) < 0){
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
  printf("minRadius:   %.2f\n", minRadius);
  printf("maxPoints:   %d\n\n", maxSize);

  Width = round2d(max.x-min.x);
  High = round2d(max.y-min.y);
  // Densidad en puntos/m^2
  Density = Npoints/(Width*High);

  // Defines a maximum number of points from which to divide the cell
  // int maxSize = (int)(minRadius*minRadius*Density);

  printf("INSERTING POINTS...\n");
  i_start = tempo_t::now();

  // Qtree qtreeIn = createQtreeF(NULL, center, maxRadius);

  // for(int i = 0; i < Npoints; i++){
  //   // insertPointF(&point_cloud[i], qtreeIn, minRadius);
  //   insertPointF2(&point_cloud[i], qtreeIn, maxSize);
  // }

  Qtree qtreeIn = parallel_qtree_creation(max_level, center, maxRadius, point_cloud, maxSize);
  // deleteQtree(parallel_qtree);
  // delete(parallel_qtree);

  printf("  CPU Reserved nodes: %zu\n", cpu_tree_nodes.load());

  i_end = tempo_t::now();
  std::cout << "INSERT CPU time elapsed: " << cast_t(i_end - i_start).count() << "ms\n";

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

  // std::vector<Lstats> vhistogram;
  std::vector<std::pair<int,int>> vhistogram;
  makeHistogram(vhistogram, qtreeIn);
  std::string out_histogram = outputTXT + ".csv";
  save_histogram(out_histogram, vhistogram, minRadius, cpu_tree_nodes, Density);
  vhistogram.clear();

  printf("Pinta la búsqueda\n");
  std::vector<int> lcpu;
  std::vector<int> lgpu;
  stageRateStats(Wsize, Overlap, Crow, Ccol, minNumPoints, lcpu, lgpu, qtreeIn, min, rate);
  printf("GPU: ");
  makeHistogram(vhistogram, lgpu);
  save_leafs(out_histogram, vhistogram);
  vhistogram.clear();
  sleep(2);
  printf("CPU: ");
  makeHistogram(vhistogram, lcpu);
  save_leafs(out_histogram, vhistogram);
  vhistogram.clear();
  lcpu.clear();
  lgpu.clear();

  printf("Pinta los niveles\n");
  std::vector<int> nlevels;
  countLevels(nlevels, qtreeIn);
  makeHistogram(vhistogram, nlevels);
  save_levels(out_histogram, vhistogram);
  vhistogram.clear();
  nlevels.clear();
  

  free(point_cloud, q);
  // free(point_cloud);
  point_cloud = NULL;

  deleteQtree(qtreeIn);
  delete(qtreeIn);

  return 0;

}
