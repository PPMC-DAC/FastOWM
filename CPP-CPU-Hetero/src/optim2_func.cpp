
Qtree_t::Qtree_t(Vector2D c, float r) : center(c), radius(r) {
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

int isLeaf(Qtree quad)
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

 void* reallocWrap(void *ptr, size_t size)
 {
     void *tmp = realloc(ptr, size);
     if(tmp == NULL)
     {
         fprintf(stderr, "Error in realloc() of size %zu\n", size);
         exit(EXIT_FAILURE);
     }
     else
     {
         ptr = tmp;
     }
     return ptr;
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

        quad->quadrants[i] = new Qtree_t(newCenter, newRadius);
    }
}

// Find the child corresponding a given point
int quadrantIdx(LpointID point, Lpoint * cloud, Qtree qtree)
{
    int child = 0;

    if(cloud[point].x >= qtree->center.x) child |= 2;
    if(cloud[point].y >= qtree->center.y) child |= 1;

    return child;
}


void insertPoint(LpointID point, Lpoint * cloud, Qtree qtree, float minRadius)
{
//    int idx = 0;

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
void tree_phase1(Qtree qtree, int maxlevel, int level, std::vector<Qtree>& n_work)
{
    if(level < maxlevel)
    {
      createQuadrants(qtree);
      for(int i = 0; i < 4; i++)
        tree_phase1(qtree->quadrants[i], maxlevel, level+1, n_work);
    }
    else
    { 
      n_work.push_back(qtree);
    }
}

/* Accumulates the points on a certain level of the tree */
void insertPointConcurrent(LpointID point, Lpoint* cloud, Qtree qtree)
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
void fillQuadrants(Lpoint* cloud, Qtree qtree, float minRadius)
{
    for(LpointID p : qtree->concurrent_points)
    {
      int idx = quadrantIdx(p, cloud, qtree);
      insertPoint(p, cloud, qtree->quadrants[idx], minRadius);
    }
    qtree->concurrent_points.clear();
}

/* Re-inserts the points when it exceeds maxNumber   */
void fillQuadrants(Lpoint* cloud, Qtree qtree, int maxNumber)
{
  for(LpointID p : qtree->concurrent_points)
    {
      int idx = quadrantIdx(p, cloud, qtree);
      insertPointMaxNum(p, cloud, qtree->quadrants[idx], maxNumber);
    }
    qtree->concurrent_points.clear();
}

//This function can be called with node_delimiter as MinRadius or MaxNumber
template<typename tree_policy> // int for MaxNumber or float for MinRadius
Qtree parallel_qtree( int level, Vector2D center, float radius, Lpoint* cloud, int Npoints, tree_policy policy )
{

  Qtree root = new Qtree_t{center, radius};

  std::vector<Qtree> n_work;
//create the "transitory" leaves up to tree-level "level" (sequential) and store these leaves in the n_work vector
  tree_phase1( root, level, 0, n_work );

  // std::cout << "  N tasks: " << n_work.size() << std::endl << std::flush;
//traverse the LiDAR points in parallel in the transitory leaves
  tbb:: parallel_for( tbb::blocked_range<int>{1, Npoints},
                      [root,cloud](tbb::blocked_range<int> r ) {
    for(int i = r.begin(); i < r.end(); i++) {
      insertPointConcurrent(i, cloud, root);
    }
  });

//Finally, traverse in parallel the transitory leaves and finish up the tree
//that hangs from them storing the points in the final leaves
  tbb::parallel_for( 0, static_cast<int>(n_work.size()), 1,
    [&](int id){

    Qtree qt = n_work[id];
    createQuadrants(qt);
    //Depending on the type of node delimiter it will call to the MinRadius or MaxNumber version
    fillQuadrants(cloud, qt, policy);
  });

  return root;
}

// Make a box with center the point and the specified radius
void makeBox(Vector2D &point, float radius, Vector2D &min, Vector2D &max)
{
    min.x = point.x - radius;
    min.y = point.y - radius;
    max.x = point.x + radius;
    max.y = point.y + radius;
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

void deleteQtree(Qtree qtree)
{
    int i;
    if(isLeaf(qtree))
    {
      qtree->points.clear();
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

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Lpoint* cloud, Qtree qtreeIn, Qtree grid, Vector2D min){

    unsigned int addMin = 0;

    #pragma omp parallel for schedule(dynamic)
    for( int jj = 0 ; jj < Ccol ; jj++ ){
       Vector2D cellCenter;
       cellCenter.y = min.y + Bsize/2 + jj*Bsize;

       for( int ii = 0 ; ii < Crow ; ii++ ){
           cellCenter.x = min.x + Bsize/2 + ii*Bsize;
           int cellPoints = 0;
           countNeighbors2D(cellCenter, cloud, grid, Bsize/2, cellPoints);
           if(cellPoints == 0){
               LpointID newmin = searchNeighborsMin(cellCenter, cloud, qtreeIn, Bsize/2, cellPoints);
               if(cellPoints>0){
                   #pragma omp critical
                   {
                      minGridIDs[addMin] = newmin;
                      addMin++;
                   }
                   // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
               }
           }
       }
    }
    return addMin;
}