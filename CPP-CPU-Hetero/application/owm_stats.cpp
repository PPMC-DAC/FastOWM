#include "../src/envi_gpu.hpp"

#include "../include/cxxopts.hpp"
#include <chrono>

int main( int argc, const char* argv[]) 
{

  // Lpoint* point_cloud = NULL; //GLOBAL VARIABLE

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
  
  cxxopts::Options options("./parallel [options]", "OWM algorithm to idetify the ground surface from LiDAR data");

  options.add_options()
          ("v,verbose", "Show data input/output", cxxopts::value<bool>()->default_value("false"))
          ("h,help", "Displays code help")
          ("i,input", "Input file name without extension where points are saved",
            cxxopts::value<std::string>()->default_value("data/INAER_2011_Alcoy"))
          ("W,Wsize", "Window size", cxxopts::value<unsigned short>()->default_value("10"))
          ("B,Bsize", "Grid size", cxxopts::value<unsigned short>()->default_value("20"))
          ("O,Overlap", "Overlap rate", cxxopts::value<double>()->default_value("0.8"))
          ("n,npes", "Number of cores", cxxopts::value<int>()->default_value("1"))
          ("l,loop", "Number of repetitions", cxxopts::value<int>()->default_value("1"))
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
  } else if ( inputTXT.find("ArzuaH.xyz") == std::string::npos &&
              inputTXT.find("AlcoyH.xyz") == std::string::npos && 
              inputTXT.find("BrionFH.xyz") == std::string::npos && 
              inputTXT.find("BrionUH.xyz") == std::string::npos ){
    printf("No header data!\n");
    exit(-1);
  }

  printf("Reading LiDAR points...\n");
  if(readXYZfile(inputTXT, point_cloud, Npoints, min, max) < 0){
    printf("Unable to read file!\n");
    exit(-1);
  }

  printf("Npoints=%lu; xmin = %.2lf; xmax = %.2lf; ymin = %.2lf; ymax = %.2lf\n",Npoints,min.x,max.x,min.y,max.y );

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


  //Qtree qtreeIn = parallel_qtree_creation(max_level, center, maxRadius, point_cloud, maxSize);
  Qtree qtreeIn= parallel_qtree_pf2(max_level, center, maxRadius, minRadius);

  i_end = tempo_t::now();
  std::cout << "INSERT CPU time elapsed: " << cast_t(i_end - i_start).count() << "ms\n";

  // Reduces the thread_private counters
  uint64_t cpu_tree_nodes = node_counter.combine([](uint64_t a, uint64_t b)
                            {return a+b;});
  printf("  CPU Reserved nodes: %zu\n", cpu_tree_nodes);


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

  // printf("Pinta la búsqueda\n");
  // std::vector<int> lcpu;
  // std::vector<int> lgpu;
  // stageRateStats(Wsize, Overlap, Crow, Ccol, minNumPoints, lcpu, lgpu, qtreeIn, min, rate);
  // printf("GPU: ");
  // makeHistogram(vhistogram, lgpu);
  // save_leafs(out_histogram, vhistogram);
  // vhistogram.clear();
  // sleep(2);
  // printf("CPU: ");
  // makeHistogram(vhistogram, lcpu);
  // save_leafs(out_histogram, vhistogram);
  // vhistogram.clear();
  // lcpu.clear();
  // lgpu.clear();

  printf("Compute the level of each leaf\n");
  std::vector<int> nlevels;
  countLevels(nlevels, qtreeIn);
  printf("Levels computed\n");
  makeHistogram(vhistogram, nlevels);
  save_levels(out_histogram, vhistogram);
  vhistogram.clear();
  nlevels.clear();
  
  freeWrap(point_cloud);

  deleteQtree(qtreeIn);

  return 0;
}
