// This version differs from optim1_qtree in the following ways:
// 1.- The function 'stage1' does not use the 'searchNeighbors2D' first to collect the points
//     inside the SW and later find the minimum if the number of points is large enough.
//     Instead, it uses findValidMin to find the minimum and the count of points inside the SW during the tree traversal.
// 2.- findValidMin does not check x,y coordinates of all points of a leaf. If the leaf's BBox is fully overlaped by the SW
//     x,y coordinates are not read. If the overlap is partial, x,y coordinate of all points have to be checked
// 3.- Stage1 in previous version used a induction variable to write the array of minimums (minIDs) so a critical section was
//     necessary for the parallel version. Now this induction is removed and minIDs can be written in parallel.
// 4.- OpenMP has been replaced by TBB in parallel stages 1 and 3
// 5.- Sorting is done with std::sort and oneapi::dpl::execution::par_unseq policy
// 6.- Tree data structure does not have id member. Now the index of the point is used as ID
// 7.- The tree is built in parallel in three phases

#include "../include/optim2_func.hpp"

int main( int argc, const char* argv[]){

  cxxopts::Options options("./owm [options]", "OWM algorithm to idetify the ground surface from LiDAR data");

  options.add_options()
          ("v,verbose", "Dump input/output data", cxxopts::value<bool>()->default_value("false"))
          ("h,help", "Help message")
          ("i,input", "Input file name without .xyz extension (LiDAR data file)",
            cxxopts::value<std::string>()->default_value("data/INAER_2011_Alcoy"))
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
  Vector2D radius = getRadius(min, max, &maxRadius); //maxRadius initialized here
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
    
  Qtree qtreeIn = parallel_qtree( level, center, maxRadius, cloud, Npoints, minRadius );
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

  // A minimum of a slide window, SW, is valid if there are at least minNumPoints inside the SW
  uint32_t minNumPoints = 0.5*Density*Wsize*Wsize;
  printf("Minium number of points in the SW to be considered: %u\n", minNumPoints);

  double Displace = round2d(Wsize*(1-Overlap));
  printf("SW x and y Displacements: %.2f\n", Displace);

  // Stage 1 parameters
  printf("\nSliding Window (SW) parameters:\n");
  ushort nCols, nRows;
  if(Overlap > 0.0) {
    nCols=(int)(round((Width+2*Wsize*Overlap)/Displace))-1;
    nRows=(int)(round((Height+2*Wsize*Overlap)/Displace))-1;
  } else {
    nCols=(int)floor(Width/Wsize)+1;
    nRows=(int)floor(Height/Wsize)+1;
  }
  printf("Number of SWs per column (number of rows): %d\n",nRows);
  printf("Number of SWs per row (number of colums):  %d\n",nCols);
  uint32_t Ncells = nCols*nRows;
  printf("Total number of OWM steps (nCols x nRows): %d\n",Ncells);

  // Stage 3 parameter
  printf("\nGrid (for stage 3) parameters:\n");
  ushort nColsg=(int)floor(Width/Bsize)+1;
  ushort nRowsg=(int)floor(Height/Bsize)+1;
  printf("Grid dimesions in %dx%d boxes: %dx%d\n\n", Bsize, Bsize, nRowsg, nColsg);
  uint32_t Ngrid = nColsg*nRowsg;
  // At most one LLP per SW and there are Ncells SWs
  int* minIDs = new int[Ncells];
  // This vector will have the points added at stage 3
  std::vector<int> minGridIDs;
  uint32_t countMin, numLLPs = 0, addMin = 0;

  tbb::tick_count t_stage, t_func;
  double t_s1, t_s2, t_s3; // execution time of stages 1, 2 and 3
  std::vector<double> results(numRuns); //To store the time results of each run
  
  for(int nR=0; nR<numRuns; nR++){

        t_func=tbb::tick_count::now();
        t_stage=tbb::tick_count::now();

        // The array minIDs stores the valid minima found 
        // The cells/SWs without a valid minimum will have a -1
        //std::fill(minIDs, minIDs+Ncells, -1);
        stage1(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, cloud, qtreeIn, min);
        t_s1=(tbb::tick_count::now()-t_stage).seconds();
        
        if(Overlap != 0){
            t_stage=tbb::tick_count::now();
            std::sort(oneapi::dpl::execution::par_unseq,minIDs,minIDs+Ncells); //std::execution::par, std::execution::par_unseq,

#ifdef DEBUG
          FILE* fileDeb;
          if((fileDeb = fopen("sortedmins.txt","w")) == nullptr){
            printf("Unable to create file!\n");
            return -1;
          }
          for(int i=0 ; i<Ncells ; i++)
              fprintf(fileDeb, "%d %.2f %.2f %.15f\n", minIDs[i], cloud[minIDs[i]].x, cloud[minIDs[i]].y,cloud[minIDs[i]].z);
          fclose(fileDeb);
#endif
            // Detect repeated ids and store them in minIDs. NumLLPs is the number of LLPs found and stored at the beggining of minIDs
            numLLPs = stage2(Ncells, minIDs);
            t_s2=(tbb::tick_count::now() - t_stage).seconds();
#ifdef DEBUG
          if((fileDeb = fopen("LLPs.txt","w")) == nullptr){
            printf("Unable to create file!\n");
            return -1;
          }
          for(int i=0 ; i<numLLPs ; i++)
              fprintf(fileDeb, "%d %.2f %.2f %.15f\n", minIDs[i], cloud[minIDs[i]].x, cloud[minIDs[i]].y,cloud[minIDs[i]].z);
          fclose(fileDeb);
#endif
        }
        else{
          countMin = std::count_if(minIDs, minIDs+Ncells, [](int i){return i>=0;} );
          numLLPs=countMin;
        }

        //printf("Number of found LLPs: %d \n\n", numLLPs);

        if(Bsize > 0){
            t_stage=tbb::tick_count::now();
            Qtree grid =  new Qtree_t( center, maxRadius) ;

            for(int i = 0; i < numLLPs; i++)
               insertPoint(minIDs[i], cloud, grid, minRadius);

            minGridIDs = stage3(Bsize, nColsg, nRowsg, cloud, qtreeIn, grid, min);
            t_s3=(tbb::tick_count::now() - t_stage).seconds();
            addMin = minGridIDs.size();
            //printf("Number of points added at stage 3: %d \n\n", addMin);
            deleteQtree(grid);
            delete(grid);
        }
        results[nR] = (tbb::tick_count::now() - t_func).seconds();
        printf("Time elapsed at STAGE 1:     %.6f s\n", t_s1);
        printf("Time elapsed at STAGE 2:     %.6f s\n", t_s2);
        printf("Time elapsed at STAGE 3:     %.6f s\n", t_s3 );
        printf("TOTAL time elapsed:     %.6f s\n", results[nR]);
        printf("Output ground seed-point cloud with %d points, %d fewer points than input cloud\n", numLLPs+addMin, Npoints-1 - numLLPs+addMin );

        printf("/////////////////////////////////////////////////////////////////////\n");
    }
    
    printf("Time elapsed at Quadtree construction:     %.6f s\n\n", time_tree);

    printf("Time of each run:  ");
    printf("  %.4lf  ", results[0]);
    for( int i=1 ; i<numRuns ; i++ ){
        printf("  %.4lf  ", results[i]);
        results[0] += results[i];
    }
    double OWMaverage=results[0]/numRuns;
    printf("\nAverage OWM: %.6lf\n", OWMaverage);
    printf("Total Tree Construction + OWM: %.6lf\n\n", OWMaverage + time_tree);
    

    for(int i = 0; i < numLLPs; i++)
       minGridIDs.push_back(minIDs[i]);
  // check results
  double correctness = check_results(goldXYZ, minGridIDs, cloud, Displace);
  if (correctness < 0) {
      printf("Unable to check results\n");
  }

  if(save_time("o2_partree.csv", inputXYZ, num_threads, minRadius, maxNumber, level,
            time_tree, OWMaverage, correctness) < 0){
    printf("Unable to create time report file!\n");
  }

#ifdef DEBUG
    // Fichero de salida
    FILE* fileMin;
    printf("Creo el fichero %s ...\n", outputXYZ);
    if((fileMin = fopen(outputXYZ.c_str(),"w")) == nullptr){
      printf("Unable to create file!\n");
      return -1;
    }

    for(int i=0 ; i<numLLPs ; i++)
      fprintf(fileMin, "%.2f %.2f %.2f\n", cloud[minIDs[i]].x, cloud[minIDs[i]].y,cloud[minIDs[i]].z);

    for(int i=0 ; i<addMin ; i++)
      fprintf(fileMin, "%.2f %.2f %.2f\n", cloud[minGridIDs[i]].x, cloud[minGridIDs[i]].y,cloud[minGridIDs[i]].z);

    //Ya no necesito mas el fichero
    if(fclose(fileMin)){
      printf("Cannot close the file\n");
      return -1;
    }
#endif
    // Free memory
    delete[] minIDs;
    minGridIDs.clear();
    freeWrap(cloud);

    deleteQtree(qtreeIn);
    delete(qtreeIn);

    return 0;
}
