#ifndef STATS_HPP

#define STATS_HPP

/* Get the density of each leaf node*/
void getDensity(std::vector<double>& vdensity, Qtree qtree) 
{
  if( isLeaf(qtree) ) { // isLeaf?
    vdensity.push_back( qtree->points.size() / (4*qtree->radius*qtree->radius) ); // density = number of points / area
  } else {
    for(int i=0; i<4; ++i) {
      getDensity(vdensity, qtree->quadrants[i]);
    }
  }
}

/* Get the radius of each leaf node*/
void getRadius(std::vector<double>& vradius, Qtree qtree) 
{
  if( isLeaf(qtree) ) { // isLeaf?
    // each time a quadrant is divided, the radius is updated
    vradius.push_back(qtree->radius);
  } else {
    for(int i=0; i<4; ++i) {
      getRadius(vradius, qtree->quadrants[i]);
    }
  }
}

/* Saves the number of points in each leaf node */
void countPoints(std::vector<int>& pCount, Qtree qtree) {

  if( isLeaf(qtree) ) { // isLeaf?
    pCount.push_back(qtree->points.size());
  } else {
    for(int i=0; i<4; ++i) {
      countPoints(pCount, qtree->quadrants[i]);
    }
  }
}


/* Saves the level of each leaf node */
void countLevels(std::vector<int>& lCount, Qtree qtree, int level = 0) 
{

  if( isLeaf(qtree) ) { // isLeaf?
    lCount.push_back(level);
  } else {
    for(int i=0; i<4; ++i) {
      countLevels(lCount, qtree->quadrants[i], level+1);
    }
  }
}


/* Creates the histogram form a tree */
void makeHistogram(std::vector<std::pair<int,int>>& histogram, Qtree qtree) 
{
  std::vector<int> pCount;
  countPoints(pCount, qtree);

  std::sort(pCount.begin(), pCount.end(), [](int& a, int& b){
          return a < b;
        });

  int idex = 0;
  int ii, jj{0}, id, reps;
  size_t endCount = pCount.size();

  for( int ii=jj ; ii<endCount ; ii=jj ) {
      //printf("ii=%d, jj=%d, endCount=%d, v[ii]=%d, v[jj]=%d\n", ii,jj,endCount, v[ii], v[jj]);
      id = pCount[ii];
      for( jj=ii+1 ; id==pCount[jj] && jj<endCount ; jj++ );
      reps = jj-ii;
      histogram.push_back({id, reps});
  }

  pCount.clear();

  // for(auto& item : histogram) {
  //   printf("%d: ", item.first);
  //   for(int i=0; i<item.second; i++) {
  //     printf("*");
  //   }
  //   printf("\n");
  // }
}


/* Creates the histogram from a vector where the count has
already been made */
template<typename data_t> // data_t is the type of the vector with the count
void makeHistogram(std::vector<std::pair<data_t,int>>& histogram, std::vector<data_t>& v) {


  std::sort(v.begin(), v.end(), [](data_t& a, data_t& b){
          return a < b;
        });

  // int idex = 0;
  data_t id;
  int ii, jj{0}, reps;
  size_t endCount = v.size();

  printf("Se han analizado %zu nodos hoja\n", endCount);

  for( int ii=jj ; ii<endCount ; ii=jj ) {
      //printf("ii=%d, jj=%d, endCount=%d, v[ii]=%d, v[jj]=%d\n", ii,jj,endCount, v[ii], v[jj]);
      id = v[ii];
      for( jj=ii+1 ; id==v[jj] && jj<endCount ; jj++ );
      reps = jj-ii;
      histogram.push_back({id, reps});
  }

}


/* Saves the leaf nodes count, with header */
template<typename policy_t> // float for MinRadius, int for MaxNumber
int save_histogram(std::string file_name, std::vector<std::pair<int,int>>& ids, policy_t policy, double density)
{
  // uint32_t size = pointList.size();
  std::fstream out(file_name, out.out | out.app);
  // std::ofstream out;
  out.precision();
  // out.open(file_name, std::ofstream::out | std::ofstream::app);
  // Point aux;
  if(!out.is_open())
    return -1;

  size_t nleaf_nodes = 0;
  for(auto& item : ids){
    nleaf_nodes += item.second;
  }

  double leaf_area;
  out << "HISTOGRAM: " << file_name << std::endl;
  out << "Number of leaf nodes: " << nleaf_nodes << std::endl;
  if constexpr(std::is_same_v<policy_t,float>){ //called with MinRadius
    out << "MinRadius: " << policy << " m" << std::endl << std::endl;
    // leaf_area = pow(2*policy, 2.);
    // out << "Area covered by MinRadius: " << leaf_area << " m^2" << std::endl << std::endl;
  }
  else{
    out << "MaxNumber: " << policy << " pts" << std::endl << std::endl;
  }

  out << "Observed: " << std::endl;

  uint32_t accumulate_n = 0;
  // double accumulate_rho = 0;
  int min_x = 999;
  int max_x = 0;
  // double min_rho = 999.;
  // double max_rho = 0.;
  for( auto& item : ids ){
    int x_i = item.first; // number of points inside this leaf-node
    /*limites*/
    if(x_i < min_x) min_x = x_i;
    if(max_x < x_i) max_x = x_i;

    int n_i = item.second; // number of leaf-nodes with this number of poinsts
    // if constexpr(std::is_same_v<policy_t,float>){ //called with MinRadius 
    //   double rho_i = x_i/leaf_area;
    //   /*limites*/
    //   if(rho_i < min_rho) min_rho = rho_i;
    //   if(max_rho < rho_i) max_rho = rho_i;
    //   accumulate_rho += rho_i*n_i;
    // }
    accumulate_n += x_i*n_i;
  }

  out << "Average number of points/leaf-node: " << accumulate_n/nleaf_nodes << std::endl;
  out << "__min: " << min_x << std::endl;
  out << "__max: " << max_x << std::endl << std::endl;
  // if constexpr(std::is_same_v<policy_t,float>){ //called with MinRadius
  //   out << "Average density: " << accumulate_rho/nleaf_nodes << std::endl;
  //   out << "__min: " << min_rho << std::endl;
  //   out << "__max: " << max_rho << std::endl << std::endl;
  // } else {
  //   out << std::endl;
  // }
  out << "Estimation (assuming equal distribution of points): " << std::endl;
  // out << "Average number of points/leaf-node: " << density*leaf_area << std::endl;
  out << "Average density: " << density << std::endl << std::endl;

  // if constexpr(std::is_same_v<policy_t,float>){ //called with MinRadius
  //   out << "x_i n_i x_i*n_i rho_i rho_i*n_i" << std::endl << std::endl;
  //   for( auto& item : ids ){
  //     int x_i = item.first; // number of points inside this leaf-node
  //     int n_i = item.second; // number of leaf-nodes with this number of poinsts
  //     double rho_i = x_i/leaf_area;
  //     out << x_i << " " << n_i << " " << x_i*n_i << " " << rho_i << " " << rho_i*n_i << std::endl;
  //   }
  // }
  // else{
    out << "x_i n_i x_i*n_i" << std::endl << std::endl;
    for( auto& item : ids ){
      int x_i = item.first; // number of points inside this leaf-node
      int n_i = item.second; // number of leaf-nodes with this number of poinsts
      out << x_i << " " << n_i << " " << x_i*n_i  << std::endl;
    }
  // }
  out << std::endl << std::endl;
  out.close();
  if(out.is_open())
    return -1;
  return 0;
}


/* Saves the leaf nodes count */
int save_leafs(std::string file_name, std::vector<std::pair<int,int>>& ids)
{
  // uint32_t size = pointList.size();
  std::fstream out(file_name, out.out | out.app);
  // std::ofstream out;
  // out.precision(std::numeric_limits<double>::digits10);
  out.precision();
  // out.open(file_name, std::ofstream::out | std::ofstream::app);
  // Point aux;
  if(!out.is_open())
    return -1;

  size_t nleaf_nodes = 0;
  for(auto& item : ids){
    nleaf_nodes += item.second;
  }
  
  out << "Number of points per leaf Histogram:"<< std::endl;

  uint64_t accumulate_n = 0;
  std::pair<int,int> max_x = {0,0};

  for( auto& item : ids ){
    int x_i = item.first; // numero de puntos en el nodo
    int n_i = item.second; // numero de veces que se repite ese número de puntos
    /*limites*/
    if(max_x.first < x_i) max_x = {x_i,n_i};

    accumulate_n += x_i*n_i;
  }

  out << "Number of leaf_nodes: " << nleaf_nodes << std::endl;
  out << "Average number of points/leaf-node: " << accumulate_n/nleaf_nodes << std::endl;
  out << "__max: " << max_x.first << " " << max_x.second << std::endl << std::endl;

  out << "x_i n_i x_i*n_i" << std::endl << std::endl;
  for( auto& item : ids ){

    int x_i = item.first; //  number of points inside this leaf-node
    int n_i = item.second; // number of leaf-nodes with this number of poinsts

    out << x_i << " " << n_i << " " << x_i*n_i << std::endl;
  }
  out << std::endl << std::endl;
  out.close();
  if(out.is_open())
    return -1;
  return 0;
}


/* Saves the radius/areas count  */
template<typename data_t>
int save_areas(std::string file_name, std::vector<std::pair<data_t,int>>& ids)
{
  // uint32_t size = pointList.size();
  std::fstream out(file_name, out.out | out.app);

  if(!out.is_open())
    return -1;

  out << "Radius HISTOGRAM:"<< std::endl;
  
  size_t nleaf_nodes = 0;
  for(auto& item : ids){
    nleaf_nodes += item.second;
  }

  out << "Number of leaf-nodes: " << nleaf_nodes << std::endl << std::endl;
  
  out << "r_i n_i a_i" << std::endl << std::endl;
  for( auto& item : ids ){

    data_t x_i = item.first; // radius
    int n_i = item.second; // number of leaf-nodes with this radius

    out << x_i << " " << n_i << " " << 4*x_i*x_i << std::endl;
  }
  out << std::endl << std::endl;
  out.close();
  if(out.is_open())
    return -1;
  return 0;
}

/* Saves the density count  */
template<typename data_t>
int save_density(std::string file_name, std::vector<std::pair<data_t,int>>& ids)
{
  // uint32_t size = pointList.size();
  std::fstream out(file_name, out.out | out.app);

  if(!out.is_open())
    return -1;

  out << "Density HISTOGRAM:"<< std::endl;
  
  size_t nleaf_nodes = 0;
  for(auto& item : ids){
    nleaf_nodes += item.second;
  }

  out << "Number of leaf-nodes: " << nleaf_nodes << std::endl << std::endl;
  
  out << "d_i n_i" << std::endl << std::endl;
  for( auto& item : ids ){

    data_t x_i = item.first; // density
    int n_i = item.second; // number of leaf-nodes with this density

    out << x_i << " " << n_i << std::endl;
  }
  out << std::endl << std::endl;
  out.close();
  if(out.is_open())
    return -1;
  return 0;
}

/* Saves the level count  */
int save_levels(std::string file_name, std::vector<std::pair<int,int>>& ids)
{
  // uint32_t size = pointList.size();
  std::fstream out(file_name, out.out | out.app);

  if(!out.is_open())
    return -1;

  out << "Levels HISTOGRAM:"<< std::endl;
  
  size_t nleaf_nodes = 0;
  for(auto& item : ids){
    nleaf_nodes += item.second;
  }

  out << "Number of leaf-nodes: " << nleaf_nodes << std::endl << std::endl;

  out << "level_i n_i" << std::endl << std::endl;
  for( auto& item : ids ){

    int x_i = item.first; // level
    int n_i = item.second; // number of leaf-nodes in this level

    out << x_i << " " << n_i << std::endl;
  }
  out << std::endl << std::endl;
  out.close();
  if(out.is_open())
    return -1;
  return 0;
}


#endif 