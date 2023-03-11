
Qtree_t::Qtree_t(Vector2D c, float r) : center(c), radius(r) {
  for(int i = 0; i < 4; i++)
    quadrants[i] = NULL;
}

QtreeConc_t::QtreeConc_t(Vector2D c, float r) : center(c), radius(r) {
  for(int i = 0; i < 4; i++)
    quadrants[i] = NULL;
}

double round2d(double z){
  return round(z*100.0)/100.0;
}

/** Calculate the radius in each axis and save the max radius of the bounding box */
Vector2D getRadius(Vector2D &min, Vector2D &max, float *maxRadius)
{
    Vector2D radii;

    radii.x = (max.x - min.x) / 2.0;
    radii.y = (max.y - min.y) / 2.0;
    if(radii.x >= radii.y)
    {
        *maxRadius = radii.x;
    }
    else
    {
        *maxRadius = radii.y;
    }

    return radii;
}

/** Calculate the center of the bounding box */
Vector2D getCenter(Vector2D &min, Vector2D &radius)
{
    Vector2D center;

    center.x = min.x + radius.x;
    center.y = min.y + radius.y;
    return center;
}

template<typename Tree>
int isLeaf(Tree quad)
{
    return quad->quadrants[0] == NULL;
}

 int isEmpty(Qtree quad)
 {
     return quad->points.size() == 0;
 }


/** Wrapper to handle malloc() nicely */
 void* mallocWrap(size_t size)
 {
     void *ptr = malloc(size);
     if(!ptr)
     {
         fprintf(stderr, "Not enough memory\n");
         exit(EXIT_FAILURE);
     }
     return ptr;
 }

template<typename pointer_t>
void freeWrap(pointer_t& ptr)
{
  free(ptr);
  ptr = nullptr;
}

void createQuadrants(Qtree quad)
{
    int i = 0;
    Vector2D newCenter;
    float newRadius = quad->radius * 0.5;

    for(i = 0; i < 4; i++)
    {
        newCenter = quad->center;
        newCenter.x += quad->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += quad->radius * (i&1 ? 0.5f : -0.5f);

        quad->quadrants[i] = new Qtree_t{newCenter, newRadius};
    }
}
void createQuadrants(QtreeConc quad)
{
    int i = 0;
    Vector2D newCenter;
    float newRadius = quad->radius * 0.5;

    for(i = 0; i < 4; i++)
    {
        newCenter = quad->center;
        newCenter.x += quad->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += quad->radius * (i&1 ? 0.5f : -0.5f);

        quad->quadrants[i] = new QtreeConc_t{newCenter, newRadius};
    }
}

// Find the child corresponding a given point
template<typename Tree>
int quadrantIdx(LpointID point, Lpoint * cloud, Tree qtree)
{
    int child = 0;

    if(cloud[point].x >= qtree->center.x) child |= 2;
    if(cloud[point].y >= qtree->center.y) child |= 1;

    return child;
}


void insertPoint(LpointID point, Lpoint * cloud, Qtree qtree, float minRadius)
{
    if(isLeaf(qtree))
    {
        if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createQuadrants(qtree);
          int idx = quadrantIdx(point, cloud, qtree);
          insertPoint(point, cloud, qtree->quadrants[idx], minRadius);
        } else {
          qtree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      int idx = quadrantIdx(point, cloud, qtree);
      insertPoint(point, cloud, qtree->quadrants[idx], minRadius);
    }
}

/* Creates nodes up to "max level" and from this level it starts
saving the pointer to the nodes in a vector that later on will be processed in parallel */
void tree_phase1(QtreeConc qtreeC, Qtree qtree , int maxlevel, int level, std::vector<std::pair<QtreeConc,Qtree>>& n_work)
{
    if(level < maxlevel)
    {
      createQuadrants(qtreeC);
      createQuadrants(qtree);
      for(int i = 0; i < 4; i++)
        tree_phase1(qtreeC->quadrants[i], qtree->quadrants[i], maxlevel, level+1, n_work);
    }
    else
    { 
      n_work.emplace_back(std::make_pair(qtreeC,qtree));
    }
}

/* Accumulates the points on a certain level of the tree */
void insertPointConcurrent(LpointID point, Lpoint* cloud, QtreeConc qtree)
{
    if(isLeaf(qtree))
    {
      // std::cout << "  Level: " << level << std::endl << std::flush;
      qtree->concurrent_points.push_back(point);
    }
    else                                // No leaf -> search the correct one
    {
      int idx = quadrantIdx(point, cloud, qtree);
      insertPointConcurrent(point, cloud, qtree->quadrants[idx]);
    }
}

void fillQuadrants(Lpoint* cloud, Qtree qtree, int maxNumber);
/* Inserts a point in the qtree creating the appropiate childs
 by setting a maximum number of points*/
void insertPointMaxNum(LpointID point, Lpoint* cloud, Qtree qtree, int maxNumber)
{
    if(isLeaf(qtree))
    {
      if(qtree->points.size() < maxNumber)
      {
        qtree->points.push_back(point);
      }
      else
      {
        createQuadrants(qtree);
        fillQuadrants(cloud, qtree, maxNumber);
        int idx = quadrantIdx(point, cloud, qtree);
        insertPointMaxNum(point, cloud, qtree->quadrants[idx], maxNumber);
      }
    }
    else                                // No leaf -> search the correct one
    {
      int idx = quadrantIdx(point, cloud, qtree);
      insertPointMaxNum(point, cloud, qtree->quadrants[idx], maxNumber);
    }
}

/* Inserts the points according to minRadius policy  */
void fillQuadrants(Lpoint* cloud, Qtree qtree, QtreeConc qtreeC, float minRadius)
{
    for(LpointID p : qtreeC->concurrent_points)
    {
      int idx = quadrantIdx(p, cloud, qtree);
      insertPoint(p, cloud, qtree->quadrants[idx], minRadius);
    }
    //qtree->concurrent_points.clear();
}

/* Re-inserts the points when it exceeds maxNumber   */
void fillQuadrants(Lpoint* cloud, Qtree qtree, QtreeConc qtreeC, int maxNumber)
{
  for(LpointID p : qtreeC->concurrent_points)
    {
      int idx = quadrantIdx(p, cloud, qtree);
      insertPointMaxNum(p, cloud, qtree->quadrants[idx], maxNumber);
    }
    //qtreeC->concurrent_points.clear();
}

void fillQuadrants(Lpoint* cloud, Qtree qtree, int maxNumber)
{
  for(LpointID p : qtree->points)
    {
      int idx = quadrantIdx(p, cloud, qtree);
      insertPointMaxNum(p, cloud, qtree->quadrants[idx], maxNumber);
    }
    //qtreeC->concurrent_points.clear();
}

//This function can be called with node_delimiter as MinRadius or MaxNumber
template<typename tree_policy> // int for MaxNumber or float for MinRadius
Qtree parallel_qtree( int level, Vector2D center, float radius, Lpoint* cloud, int Npoints, tree_policy policy )
{

  QtreeConc rootConc = new QtreeConc_t{center, radius};
  Qtree root = new Qtree_t{center, radius};

  std::vector<std::pair<QtreeConc,Qtree>> n_work;
//create the "transitory" leaves up to tree-level "level" (sequential) and store these leaves in the n_work vector
  tree_phase1( rootConc, root, level, 0, n_work );

  // std::cout << "  N tasks: " << n_work.size() << std::endl << std::flush;
//traverse the LiDAR points in parallel in the transitory leaves
  tbb:: parallel_for( tbb::blocked_range<int>{1, Npoints},
                      [rootConc,cloud](tbb::blocked_range<int> r ) {
    for(int i = r.begin(); i < r.end(); i++) {
      insertPointConcurrent(i, cloud, rootConc);
    }
  });

//Finally, traverse in parallel the transitory leaves and finish up the tree
//that hangs from them storing the points in the final leaves
  tbb::parallel_for( 0, static_cast<int>(n_work.size()), 1,
    [&](int id){

    Qtree qt = n_work[id].second;
    QtreeConc qtc = n_work[id].first;
    createQuadrants(qt);
    //Depending on the type of node delimiter it will call to the MinRadius or MaxNumber version
    fillQuadrants(cloud, qt, qtc, policy);
  });
  deleteQtree(rootConc);
  delete(rootConc);
  return root;
}

// Makes a square box centered at center with the specified radius
void makeBox(Vector2D &point, float radius, Vector2D &min, Vector2D &max)
{
    min.x = point.x - radius;
    min.y = point.y - radius;
    max.x = point.x + radius;
    max.y = point.y + radius;
}

/* Makes a box centered at point center (can be rectangular if radiusX != radiusY) */
void makeBox(Vector2D& center, double radiusX, double radiusY, Vector2D& min, Vector2D& max)
{
    min.x = center.x - radiusX;
    min.y = center.y - radiusY;
    max.x = center.x + radiusX;
    max.y = center.y + radiusY;
}

int insideBox2D(Lpoint* cloud, LpointID point, Vector2D &min, Vector2D &max)
{
    double x=cloud[point].x;
    double y=cloud[point].y;
    if(x > min.x && y > min.y && x < max.x && y < max.y) return 1;
    return 0;
}

int boxInside2D(Vector2D &boxMin, Vector2D &boxMax, Qtree quad)
{
    if(quad->center.x + quad->radius > boxMax.x ||
       quad->center.y + quad->radius > boxMax.y)
        return 0;

    if(quad->center.x - quad->radius < boxMin.x ||
       quad->center.y - quad->radius < boxMin.y)
        return 0;

    return 1;
}

int boxOverlap2D(Vector2D &boxMin, Vector2D &boxMax, Qtree quad)
{
    if(quad->center.x + quad->radius < boxMin.x ||
       quad->center.y + quad->radius < boxMin.y)
        return 0;

    if(quad->center.x - quad->radius > boxMax.x ||
       quad->center.y - quad->radius > boxMax.y)
        return 0;

    return 1;
}

int boxTotalOverlap2D(Vector2D &boxMin, Vector2D &boxMax, Qtree &quad)
{
    if(quad->center.x + quad->radius < boxMax.x ||
       quad->center.y + quad->radius < boxMax.y)
        return 0;

    if(quad->center.x - quad->radius > boxMin.x ||
       quad->center.y - quad->radius > boxMin.y)
        return 0;

    return 1;
}

template<typename Tree>
void deleteQtree(Tree qtree)
{
    int i;
    if(isLeaf(qtree))
    {
      if constexpr (std::is_same_v<Tree,Qtree>) qtree->points.clear();
      else qtree->concurrent_points.clear();
    } else {
        for(i = 0; i < 4; i++) {
            // Check
            deleteQtree(qtree->quadrants[i]);
            delete(qtree->quadrants[i]);
        }
    }
}

void findValidMin(Lpoint* cloud, Qtree qtree, Vector2D &boxMin, Vector2D &boxMax, int &numInside, LpointID &minidx)
{
    if(isLeaf(qtree))
    {
      if(boxInside2D(boxMin, boxMax, qtree)){
        //if(minidx==-1 && qtree->points.size()) minidx=qtree->points[0];
        for(LpointID p : qtree->points) {
          if (cloud[p].z < cloud[minidx].z) {
              minidx = p;
          }
        }
        numInside+=qtree->points.size(); //passed by reference
      } else {
        for(LpointID p : qtree->points) { 
            double x=cloud[p].x;
            double y=cloud[p].y;
            if(x > boxMin.x && y > boxMin.y && x < boxMax.x && y < boxMax.y)
            { 
              if (cloud[p].z < cloud[minidx].z) minidx = p;
              numInside++;
            }
        }
      }

    } else {
        for(int i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(boxMin, boxMax, qtree->quadrants[i]))
                continue;
            else {
                findValidMin(cloud, qtree->quadrants[i], boxMin, boxMax, numInside, minidx);
            }
        }
    }
}

LpointID searchNeighborsMin(Vector2D &SW_center, Lpoint* cloud, Qtree qtree, float radius, int &numInside)
{
    Vector2D boxMin, boxMax;
    LpointID minidx = 0; // In position 0 we have this point: Lpoint{0.0, 0.0, std::numeric_limits<double>::max()};
    numInside = 0;
    makeBox(SW_center, radius, boxMin, boxMax); //updates boxMin,boxMax to be de BBox of SW with center and radius

    findValidMin(cloud, qtree, boxMin, boxMax, numInside, minidx);
    return minidx; 
}

LpointID searchNeighborsMin(Vector2D &SW_center, Lpoint* cloud, Qtree qtree, float radiusX, float radiusY, int &numInside)
{
    Vector2D boxMin, boxMax;
    LpointID minidx = 0; // In position 0 we have this point: Lpoint{0.0, 0.0, std::numeric_limits<double>::max()};
    numInside = 0;
    makeBox(SW_center, radiusX, radiusY, boxMin, boxMax); //updates boxMin,boxMax to be de BBox of SW with center and radius

    findValidMin(cloud, qtree, boxMin, boxMax, numInside, minidx);
    return minidx; 
}

void countNeighbors(Lpoint* cloud, Qtree qtree, Vector2D &boxMin, Vector2D &boxMax, int &numInside)
{
    int i;

    if(isLeaf(qtree))
    {
        if(boxInside2D(boxMin, boxMax, qtree)){
            numInside+=qtree->points.size();
        }
        else{
          for(LpointID p : qtree->points) {
            double x=cloud[p].x;
            double y=cloud[p].y;
            if(x > boxMin.x && y > boxMin.y && x < boxMax.x && y < boxMax.y) numInside++;
          }
        }
    } else {
        for(i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(boxMin, boxMax, qtree->quadrants[i]))
              continue;
            else {
              countNeighbors(cloud, qtree->quadrants[i], boxMin, boxMax, numInside);
            }
        }
    }
}

void countNeighbors2D(Vector2D &point, Lpoint* cloud, Qtree qtree, float radius, int &numInside)
{
    Vector2D boxMin, boxMax;

    numInside = 0;
    makeBox(point, radius, boxMin, boxMax);

    countNeighbors(cloud, qtree, boxMin, boxMax, numInside);
}

void stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Lpoint* cloud, Qtree qtreeIn, Vector2D min){

  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  tbb::parallel_for(tbb::blocked_range<int>{0,Crow*Ccol},
                      [&](tbb::blocked_range<int> r ) {
        for( int step = r.begin() ; step < r.end() ; step++ ){
                int ii=step/Ccol, jj=step%Ccol;
                Vector2D cellCenter={initX + ii*Displace, initY + jj*Displace};
                int cellPoints = 0;
                LpointID newmin = searchNeighborsMin(cellCenter, cloud, qtreeIn, Wsize/2, cellPoints);
                //printf("Step: %d.%d; Min id: %.2f; cellPoints: %d\n",ii,jj,newmin, cellPoints);
                if(cellPoints >= minNumPoints ) minIDs[step] = newmin;
                else minIDs[step] = -1;
        }
  });
}

void stage1memoA(ushort Wsize, double Overlap, ushort nCols, ushort nRows,
  ushort minNumPoints, int* minIDs, Lpoint* cloud, Qtree qtreeIn, Vector2D min){

  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  tbb::parallel_for(tbb::blocked_range2d<int,int>{0,nRows,0,nCols},
                      [&](tbb::blocked_range2d<int,int> r ) {
      LpointID newmin = 0;
      Vector2D cellCenter;
      Vector2D boxMax, boxMin;
      int je = r.rows().end();
      int ie = r.cols().end();
      for( int jj = r.rows().begin(); jj < je; ++jj){
        int cellPoints = 0;
        cellCenter.y = initY + jj*Displace;
        int jstart = jj*nCols;
        for(int ii = r.cols().begin() ; ii < ie ; ii++ ){
          cellCenter.x=initX + ii*Displace;
          makeBox(cellCenter, Wsize*0.5, boxMin, boxMax);
          //Memoization: if the previous minimum is inside the box, we don't need to search in the whole SW
          if(insideBox2D(cloud, newmin, boxMin, boxMax)){
            //Smaller area in which we have to search now (only the one not visited in previous step)
            Vector2D smallerCellCenter = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};
            int old_cellPoints = cellPoints; //remember the points counted in the previous step
            // Search in the smaller area Displace x Wsize centered at new cellCenter
            LpointID tmpmin = searchNeighborsMin(smallerCellCenter, cloud, qtreeIn, Displace*0.5, Wsize*0.5, cellPoints);
            // We're assuming homogeneous densities around the SW to compute the number of points inside the SW.
            /*If we need more precission, we could count with countMin() the number of points inside SW*/
            cellPoints += static_cast<int>(old_cellPoints * Overlap);
            // If the min of the reduced area is smaller, use that one
            if(cloud[tmpmin].z < cloud[newmin].z) newmin = tmpmin;
          } else {
            //Otherwise we have to search in the whole SW of Wsize x Wsize
            newmin = searchNeighborsMin(cellCenter, cloud, qtreeIn, Wsize/2, Wsize/2, cellPoints);
          }
          //printf("Step: %d; Min id: %d; cellPoints: %d\n",step, newmin, cellPoints);
          if(cellPoints >= minNumPoints) minIDs[jstart + ii] = newmin;
          else minIDs[jstart + ii] = -1;
        }
      }
  });
}


//Receives a sorted list of minIDs with -1s at the beginning due to the SWs that do not have an LLP
unsigned int stage2(unsigned int Ncells, int* minIDs){

  int i=0;
  while(minIDs[i]<0 && i<Ncells) i++; //Skip -1s
  unsigned int counter = 0;
  int jj=0;
  for( int ii=i ; ii<Ncells ; ii=jj ){
      int id = minIDs[ii];
      for( jj=ii+1 ; id==minIDs[jj] && jj<Ncells ; jj++ );
      if(jj-ii > 1){
          minIDs[counter]=id;
          counter++;
      }
  }
  return counter;
}

std::vector<int> stage3(unsigned short Bsize, unsigned short nCols, unsigned short nRows,
           Lpoint* cloud, Qtree qtreeIn, Qtree grid, Vector2D min){

    double initX =  min.x + Bsize/2;
    double initY = min.y + Bsize/2;

    std::vector<int> init;
    return tbb::parallel_reduce(
      tbb::blocked_range2d<int,int>{0,nRows,0,nCols},
      init,
      [&](const tbb::blocked_range2d<int,int>& r, std::vector<int> v) {

        Lpoint newmin;
        Vector2D cellCenter;
        int cellPoints;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {
            cellCenter.y = initY + jj*Bsize;
            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){
                cellCenter.x = initX + ii*Bsize;
                countNeighbors2D(cellCenter, cloud, grid, Bsize/2, cellPoints);
                if(cellPoints == 0){
                    LpointID newmin = searchNeighborsMin(cellCenter, cloud, qtreeIn, Bsize/2, cellPoints);
                    if(cellPoints>0) v.push_back(newmin);
                }
           }
       }
       return v;
    },
    [](std::vector<int> a, std::vector<int> b) {
        a.insert(a.end(), b.begin(), b.end());
        return a;
    });
}

/* Reads the point cloud in .xyz format */
int readXYZfile(std::string filename, Lpoint* & point_cloud, uint32_t & Npoints, Vector2D &min, Vector2D &max) 
{
  //Long file names correspond to xyz files without header data, so we provide the info here:
  if( filename.find("INAER_2011_Alcoy.xyz") != std::string::npos ){ // Alcoy mini
    Npoints = 2772832; min={715244.96,4286623.63}; max={716057.75,4287447.70};
  } else if( filename.find("INAER_2011_Alcoy_Core.xyz") != std::string::npos ){ // Alcoy
    Npoints = 20380212; min={714947.98, 4286501.93}; max={716361.06, 4288406.23};
  } else if( filename.find("BABCOCK_2017_Arzua_3B.xyz") != std::string::npos ){ //Arzua
    Npoints = 40706503; min={568000.00, 4752320.00}; max={568999.99, 4753319.99};
  } else if( filename.find("V21_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion forestal
    Npoints = 42384876; min={526964.093, 4742610.292}; max={527664.647, 4743115.738};
  } else if( filename.find("V19_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion urban
    Npoints = 48024480; min={526955.908, 4742586.025}; max={527686.445, 4743124.373};
  } else if( filename.find("sample24.xyz") != std::string::npos ){
    Npoints = 7492; min={513748.12, 5403124.76}; max={513869.97, 5403197.20};
    //These files with short name *H.xyz DO HAVE header info inside the file:
  } else if ( filename.find("ArzuaH.xyz") == std::string::npos &&
              filename.find("AlcoyH.xyz") == std::string::npos && 
              filename.find("BrionFH.xyz") == std::string::npos && 
              filename.find("BrionUH.xyz") == std::string::npos ){
    printf("No header data!\n");
    exit(-1);
  }
  FILE* fileXYZ;
  if((fileXYZ = fopen(filename.c_str(),"r")) == NULL){
    printf("Unable to open file!\n");
    return -1;
  }
  if ( filename.find("ArzuaH.xyz") != std::string::npos || filename.find("AlcoyH.xyz") != std::string::npos || 
       filename.find("BrionFH.xyz") != std::string::npos || filename.find("BrionUH.xyz") != std::string::npos ){
    printf("Read header...\n");
    if(fscanf(fileXYZ, "%d\n%lf\n%lf\n%lf\n%lf\n",&Npoints, &min.x, &max.x, &min.y, &max.y) < 5){
        printf("Imposible to read header values\n");
        return -1;
    }
  }
  printf("Npoints=%d; xmin = %.2lf; xmax = %.2lf; ymin = %.2lf; ymax = %.2lf\n",Npoints,min.x,max.x,min.y,max.y );
  Npoints++; // we insert at position 0 an artificial point to compare: {0.0, 0.0, std::numeric_limits<double>::max()}

  // Allocate memory for the LiDAR points
  point_cloud = static_cast<Lpoint*>(mallocWrap(Npoints * sizeof(Lpoint)));

  printf("Reading LiDAR points...\n");
  point_cloud[0]={0.0, 0.0, std::numeric_limits<double>::max()};
  for(int i=1; i<Npoints ; i++){
    if(fscanf(fileXYZ, "%lf %lf %lf",&point_cloud[i].x,&point_cloud[i].y,&point_cloud[i].z) < 3){
      printf("Error reading values\n");
      return -1;
    }
    while(fgetc(fileXYZ)!='\n');
  }

  if(fclose(fileXYZ)){
    printf("Cannot close the file\n");
    return -1;
  }
  return 0;
}

/* Checks results by comparing with GOLD result */
double check_results(std::string filename, std::vector<int>& ids, Lpoint* cloud, float Displace)
{
  std::string line;
  // char delim;
  Vector3D aux;
  uint32_t npoints = 0;
  std::vector<Vector3D> vgold, v;

  std::ifstream input(filename);
  if (input.fail()) {
      std::cout << "The GOLD point list doesn't exits" << std::endl;
      return -1.;
  }

  while(getline(input,line)) {
    std::stringstream s(line);
    // s >> aux.x >> delim >> aux.y >> delim >> aux.z;
    s >> aux.x >> aux.y >> aux.z;
    vgold.push_back(aux);
    npoints++;
  }

  (input).close();
  if((input).is_open())
    return -1.;

  std::sort(oneapi::dpl::execution::par_unseq, vgold.begin(), vgold.end(), [](const Vector3D a, const Vector3D b){
    return a.z < b.z;
  });

  // for(auto p1 : v){
  //   printf("%lf %lf %lf\n", p1.x, p1.y, p1.z);
  // }

  for( int i : ids ){
    aux.x = cloud[i].x;
    aux.y = cloud[i].y;
    aux.z = cloud[i].z;
    v.push_back(aux);
  }

  std::sort(oneapi::dpl::execution::par_unseq, v.begin(), v.end(), [](const Vector3D a, const Vector3D b){
    return a.z < b.z;
  });

  int count = 0;
  // float Displace = 1.0;
  Vector2D boxMin, boxMax;
  for(auto p : v){
    Vector2D aux = {p.x, p.y};
    makeBox(aux, Displace*0.5, boxMin, boxMax);
    for(auto pg : vgold){
      if(pg.x > boxMin.x && pg.y > boxMin.y && pg.x < boxMax.x && pg.y < boxMax.y){
        if(fabs(p.z - pg.z) < 0.01){
          count++;
          break;
          }
        }
    }
  }

  double ratio = count/((double)(npoints))*100.0;
  printf("%d points correct; %.2f%%\n", count, ratio);

  return ratio;

  // return std::equal(v.begin(), v.end(), v2.begin(), [](const Vector3D a, const Vector3D b){
  //   return fabs(a.z - b.z) < 0.01;
  // });
}

/* Saves the execution time and other metrics of the algorithm  */
int save_time(const char* exec_name, std::string file_name, std::string map_name, int numthreads, 
  float minRadius, int maxNumber, int level, double tree_time,  
  double owm_time, double correctness)
{
  bool newFile=false;
  if(!std::filesystem::exists(file_name)) newFile=true;

  // uint32_t size = pointList.size();
  std::fstream out(file_name, out.out | out.app);
  // std::ofstream out;
  // out.precision(std::numeric_limits<double>::digits10);
  // out.open(file_name);
  // Point aux;
  if(!out.is_open())
    return -1;

  if(newFile){
    out << "ExecName InputFile NumThreads MinRadius MaxNumber Level TreeTime OwmTime Correctness\n";
  }
  out << exec_name << " " << map_name << " " << numthreads << " ";
  // out.precision(3);
  out << std::defaultfloat << minRadius << " " << maxNumber << " " << level << " ";
  // out.precision(6);
  out << std::fixed << tree_time << " " <<  owm_time << " ";
  // out.precision(1);
  out << std::defaultfloat  << correctness << std::endl;

  out.close();
  if(out.is_open())
    return -1;
  return 0;
}