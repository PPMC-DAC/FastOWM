
#include "../include/optim3_func.hpp"
#include "../include/stats_func.hpp"

int main( int argc, const char* argv[]){

  cxxopts::Options options("./owm [options]", "OWM algorithm to idetify the ground surface from LiDAR data");

  options.add_options()
          ("v,verbose", "Dump input/output data", cxxopts::value<bool>()->default_value("false"))
          ("h,help", "Help message")
          ("i,input", "Input file name without .xyz extension (LiDAR data file)",
            cxxopts::value<std::string>()->default_value("data/INAER_2011_Alcoy"))
          ("results", "Output results file",
            cxxopts::value<std::string>()->default_value("tbb_histogram.csv"))
          ("W,Wsize", "Sliding Window size", cxxopts::value<uint32_t>()->default_value("10"))
          ("B,Bsize", "Grid size", cxxopts::value<uint32_t>()->default_value("20"))
          ("O,Overlap", "Overlap ratio", cxxopts::value<double>()->default_value("0.80"))
          ("n,num_threads", "Number of threads", cxxopts::value<int>()->default_value("1"))
          ("l,loop", "Number of runs of the OWM algorithm", cxxopts::value<int>()->default_value("1"))
          ("r,radius", "Value of minRadius (tree cutoff policy)", cxxopts::value<float>()->default_value("0.1"))
          ("s,size", "Value of maxNumber (max number of points per leaf-node)", cxxopts::value<int>()->default_value("32"))
          ("L,level", "Tree level at which the tree creation becomes parallel", 
            cxxopts::value<int>()->default_value("5"));

  auto parameters = options.parse(argc, argv);

  uint32_t Wsize = parameters["Wsize"].as<uint32_t>();
  uint32_t Bsize = parameters["Bsize"].as<uint32_t>();
  double Overlap = parameters["Overlap"].as<double>();
  size_t num_threads = parameters["num_threads"].as<int>();
  int numRuns = parameters["loop"].as<int>();
  float minRadius = parameters["radius"].as<float>();
  int maxNumber = parameters["size"].as<int>();
  int level = parameters["level"].as<int>();

  std::string inputXYZ = parameters["input"].as<std::string>() + ".xyz";
  std::string outputXYZ = parameters["input"].as<std::string>() + "_salida.xyz";
  std::string goldXYZ = parameters["input"].as<std::string>() + "_salidaGOLD.xyz";

  // Results file
  std::string resultsCSV = parameters["results"].as<std::string>();

  if (parameters.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
  }

  //omp_set_num_threads(num_threads);
  tbb::global_control c{tbb::global_control::max_allowed_parallelism, num_threads};

  printf("Input.txt: %s ---> EX. CON %ld CORES\n", inputXYZ.c_str(), num_threads);
  uint32_t Npoints=0; // Number of points of the LiDAR point cloud
  Vector2D  min, max; //min and max (x,y) coordinates of each cloud

  Lpoint* cloud=nullptr; //Allocated and filled inside factory readXYZfile()
  printf("Reading LiDAR points...\n");
  if(readXYZfile(inputXYZ, cloud, Npoints, min, max) < 0){
    printf("Unable to read file!\n");
    exit(-1);
  }
  float maxRadius = 0.0; // highest distance (x or y direction) of the cloud area
  Vector2D radius = getRadius(min, max, maxRadius); //maxRadius by ref initialized here
  Vector2D center = getCenter(min, radius);
  printf("QTREE PARAMETERS:\n");
  printf("MaxRadius:  %.3f\n", maxRadius);
  printf("Center:     %.2f , %.2f\n", center.x,center.y);
  printf("Radius:     %.2f , %.2f\n", radius.x,radius.y);
  printf("BUILDING QTREE...\n");

  printf("Building quadtree with minRadius: %g\n", minRadius);
  tbb::tick_count t_qtree=tbb::tick_count::now();

  // Alternative that construct the tree sequentially: 
  // Qtree qtreeIn = new Qtree_t( center, maxRadius );
  // for(int i = 1; i <= Npoints; i++)
  //    insertPoint(i, cloud, qtreeIn, minRadius);
  #ifdef MAXNUMBER
  Qtree qtreeIn = parallel_qtree( level, center, maxRadius, cloud, Npoints, maxNumber );
  #else
  Qtree qtreeIn = parallel_qtree( level, center, maxRadius, cloud, Npoints, minRadius );
  #endif  
  storeMinAndNumPoints(cloud, qtreeIn, 0);
  double time_tree = (tbb::tick_count::now()-t_qtree).seconds();
  
  double Width = round2d(max.x-min.x);
  double Height = round2d(max.y-min.y);
  double Density = (Npoints-1)/(Width*Height);
  printf("CLOUD PARAMETERS:\n");
  printf("Number of points: %d\n",Npoints-1);
  printf("Width:   %.2lf\n",Width);
  printf("Height:  %.2lf\n",Height);
  printf("Density: %.3lf\n",Density);

  printf("\nSize of the sliding window (SW): %u\n", Wsize);
  printf("Size of the grid cell: %u\n", Bsize);
  printf("Overlap:               %.2f\n", Overlap);

  std::vector<std::pair<int,int>> vhistogram;
  makeHistogram(vhistogram, qtreeIn);
  #ifdef MAXNUMBER
    save_histogram(resultsCSV, vhistogram, maxNumber, Density);
  #else
    save_histogram(resultsCSV, vhistogram, minRadius, Density);
  #endif
  vhistogram.clear();

  printf("Get the radius of each leaf node\n");
  std::vector<double> vradius;
  getRadius(vradius, qtreeIn);
  printf("Radius computed\n");
  std::vector<std::pair<double,int>> ahistogram;
  makeHistogram(ahistogram, vradius);
  save_areas(resultsCSV, ahistogram);
  vradius.clear();
  ahistogram.clear();

  printf("Get the density of each leaf node\n");
  std::vector<double> vdensity;
  getDensity(vdensity, qtreeIn);
  printf("Density computed\n");
  std::vector<std::pair<double,int>> dhistogram;
  makeHistogram(dhistogram, vdensity);
  save_density(resultsCSV, dhistogram);
  vdensity.clear();
  dhistogram.clear();

  printf("Compute the level of each leaf\n");
  std::vector<int> nlevels;
  countLevels(nlevels, qtreeIn);
  printf("Levels computed\n");
  makeHistogram(vhistogram, nlevels);
  save_levels(resultsCSV, vhistogram);
  vhistogram.clear();
  nlevels.clear();
  
  freeWrap(cloud);

  deleteQtree(qtreeIn);

  return 0;
}
