#ifndef STATS_HPP

#define STATS_HPP

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
void makeHistogram(std::vector<std::pair<int,int>>& histogram, std::vector<int>& v) {


  std::sort(v.begin(), v.end(), [](int& a, int& b){
          return a < b;
        });

  int idex = 0;
  int ii, jj{0}, id, reps;
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
int save_histogram(std::string file_name, std::vector<std::pair<int,int>>& ids, float minRadius, uint64_t num_nodes, 
  double density)
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

  double leaf_area = pow(2*minRadius, 2.);
  
  out << "HISTOGRAM: " << file_name << std::endl;
  out << "N_puntos: " << Npoints << std::endl;
  out << "N_nodos: " << num_nodes << std::endl;
  out << "N_nodos_hoja: " << nleaf_nodes << std::endl;
  out << "MinRadius: " << minRadius << " m" << std::endl;
  out << "Area_minima: " << leaf_area << " m^2" << std::endl << std::endl;
  out << "OBSERVADO: " << std::endl;

  uint32_t accumulate_n = 0;
  double accumulate_rho = 0;
  int min_x = 999;
  int max_x = 0;
  double min_rho = 999.;
  double max_rho = 0.;
  for( auto& item : ids ){
    int x_i = item.first; // numero de puntos en el nodo
    /*limites*/
    if(x_i < min_x) min_x = x_i;
    if(max_x < x_i) max_x = x_i;

    int n_i = item.second; // numero de veces que se repite ese número de puntos

    double rho_i = x_i/leaf_area;
    /*limites*/
    if(rho_i < min_rho) min_rho = rho_i;
    if(max_rho < rho_i) max_rho = rho_i;

    accumulate_n += x_i*n_i;
    accumulate_rho += rho_i*n_i;
  }

  // printf("densidad media :::: %g", accumulate_rho/nleaf_nodes);

  out << "N_medio_pts/nodo: " << accumulate_n/nleaf_nodes << std::endl;
  out << "__min: " << min_x << std::endl;
  out << "__max: " << max_x << std::endl;
  out << "Densidad_media: " << accumulate_rho/nleaf_nodes << std::endl;
  out << "__min: " << min_rho << std::endl;
  out << "__max: " << max_rho << std::endl << std::endl;
  out << "ESTIMADO: " << std::endl;
  out << "N_medio_pts/nodo: " << density*leaf_area << std::endl;
  out << "Densidad_media: " << density << std::endl << std::endl;

  out << "x_i n_i x_i*n_i rho_i rho_i*n_i" << std::endl << std::endl;
  for( auto& item : ids ){
    int x_i = item.first; // numero de puntos en el nodo
    int n_i = item.second; // numero de veces que se repite ese número de puntos
    double rho_i = x_i/leaf_area;
    out << x_i << " " << n_i << " " << x_i*n_i << " " << rho_i << " " << rho_i*n_i << std::endl;
    // out << item.first << " " << item.second << std::endl;
    // fprintf(out, "%.2f %.2f %.2f\n", (*pointer)[i].x, (*pointer)[i].y, (*pointer)[i].z);
  }
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
  
  out << "HISTOGRAM_DE_BUSQUEDA:"<< std::endl;

  uint64_t accumulate_n = 0;
  std::pair<int,int> max_x = {0,0};

  for( auto& item : ids ){
    int x_i = item.first; // numero de puntos en el nodo
    int n_i = item.second; // numero de veces que se repite ese número de puntos
    /*limites*/
    if(max_x.first < x_i) max_x = {x_i,n_i};

    accumulate_n += x_i*n_i;
  }

  out << "Nodos_hoja_analizados: " << nleaf_nodes << std::endl;
  out << "N_medio_pts/nodo: " << accumulate_n/nleaf_nodes << std::endl;
  out << "__max: " << max_x.first << " " << max_x.second << std::endl << std::endl;

  out << "x_i n_i x_i*n_i" << std::endl << std::endl;
  for( auto& item : ids ){

    int x_i = item.first; // numero de puntos en el nodo
    int n_i = item.second; // numero de veces que se repite ese número de puntos

    out << x_i << " " << n_i << " " << x_i*n_i << std::endl;
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

  
  out << "HISTOGRAM_DE_NIVELES:"<< std::endl;
  
  size_t nleaf_nodes = 0;
  for(auto& item : ids){
    nleaf_nodes += item.second;
  }

  out << "Nodos_hoja_analizados: " << nleaf_nodes << std::endl << std::endl;

  out << "x_i n_i" << std::endl << std::endl;
  for( auto& item : ids ){

    int x_i = item.first; // numero de puntos en el nodo
    int n_i = item.second; // numero de veces que se repite ese número de puntos

    out << x_i << " " << n_i << std::endl;
  }
  out << std::endl << std::endl;
  out.close();
  if(out.is_open())
    return -1;
  return 0;
}


#endif 