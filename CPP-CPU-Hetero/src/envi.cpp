#include "../include/envi.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/concurrent_vector.h"
#include "tbb/blocked_range2d.h"
#include "tbb/task.h"

// #include "tbb/cache_aligned_allocator.h"
// #include "tbb/tbb.h"
// #include "tbb/scalable_allocator.h" // use -ltbbmalloc
// #include "tbb/tbbmalloc_proxy.h"
// #include <numeric>
#include <atomic>
// #include <mutex>

#include <fstream>
#include <sstream>


Qtree_t::Qtree_t(Vector2D c, float r) : center(c), radius(r) {
  for(int i = 0; i < 4; i++)
    quadrants[i] = NULL;
}

double round2d(double z){
  return round(z*100.0)/100.0;
}

// Make a box with center the point and the specified radius
void makeBox(Vector2D *point, float radius, Vector2D *min, Vector2D *max)
{
    // printf("Radio: %.2f\n", radius);
    // printf("Centro: [ %.2lf, %.2lf]\n",point->x,point->y );
    min->x = point->x - radius;
    min->y = point->y - radius;
    max->x = point->x + radius;
    max->y = point->y + radius;

}

void makeBox(Vector2D *point, double radiusX, double radiusY, Vector2D *min, Vector2D *max)
{

    min->x = point->x - radiusX;
    min->y = point->y - radiusY;

    max->x = point->x + radiusX;
    max->y = point->y + radiusY;

}

int boxInside2D(Vector2D boxMin, Vector2D boxMax, Qtree qt)
{
    if(qt->center.x + qt->radius > boxMax.x ||
       qt->center.y + qt->radius > boxMax.y)
        return 0;

    if(qt->center.x - qt->radius < boxMin.x ||
       qt->center.y - qt->radius < boxMin.y)
        return 0;

    return 1;
}

int insideBox2D(Vector3D *point, Vector2D min, Vector2D max)
{
    if(point->x > min.x && point->y > min.y)
    {
        if(point->x < max.x && point->y < max.y)
        {
            return 1;
        }
    }

    return 0;
}

int insideBox2D(Lpoint *point, Vector2D min, Vector2D max)
{
    if(point->x > min.x && point->y > min.y)
    {
        if(point->x < max.x && point->y < max.y)
        {
            return 1;
        }
    }

    return 0;
}


// int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, Qtree_t* qt)
int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, Qtree qt)
{
    if(qt->center.x + qt->radius < boxMin.x ||
       qt->center.y + qt->radius < boxMin.y)
        return 0;

    if(qt->center.x - qt->radius > boxMax.x ||
       qt->center.y - qt->radius > boxMax.y)
        return 0;

    return 1;
}

int boxTotalOverlap2D(Vector2D boxMin, Vector2D boxMax, Qtree qt)
{
    if(qt->center.x + qt->radius < boxMax.x ||
       qt->center.y + qt->radius < boxMax.y)
        return 0;

    if(qt->center.x - qt->radius > boxMin.x ||
       qt->center.y - qt->radius > boxMin.y)
        return 0;

    return 1;
}

int isLeaf(Qtree qt)
{
    return qt->quadrants[0] == NULL;
}

/** Calculate the radius in each axis and save the max radius of the bounding box */
Vector2D getRadius(Vector2D min, Vector2D max, float *maxRadius)
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
Vector2D getCenter(Vector2D min, Vector2D radius)
{
    Vector2D center;

    center.x = min.x + radius.x;
    center.y = min.y + radius.y;

    return center;
}

// Find the child index corresponding to a given point
int quadrantIdx(Lpoint *point, Qtree qtree)
{
    int child = 0;

    if(point->x >= qtree->center.x) child |= 2;
    if(point->y >= qtree->center.y) child |= 1;

    return child;
}


// Lpoint findValidMin(Qtree qtree, const std::vector<std::function<bool(Vector2D, float)>>& conditions)
Lpoint findValidMin(Qtree qtree, Vector2D* boxMin, Vector2D* boxMax, int* numInside)
{
    Lpoint tmp, min = {0,0.0,0.0,std::numeric_limits<double>::max()};

    if(isLeaf(qtree))
    {
      if(boxInside2D(*boxMin, *boxMax, qtree)){
        for(Lpoint* p : qtree->points) {
          if (p->z < min.z) {
              min = *p;
          }
          (*numInside)++;
        }
      } else {
        for(Lpoint* p : qtree->points) {
          if(insideBox2D(p, *boxMin, *boxMax))
          {
            if (p->z < min.z) {
                min = *p;
            }
            (*numInside)++;
          }
        }
      }
    } else {
        for(int i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, qtree->quadrants[i]))
                continue;
            else {
                tmp = findValidMin(qtree->quadrants[i], boxMin, boxMax, numInside);
                if (tmp.z < min.z) {
                    min = tmp;
                }
            }
        }
    }

    return min;
}

Lpoint searchNeighborsMin(Vector2D* point, Qtree qtree, float radius, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMin(qtree, &boxMin, &boxMax, numNeighs);
}

Lpoint findValidMin(Qtree qtree, Vector2D* boxMin, Vector2D* boxMax, int* numInside, int level)
{
    // Lpoint tmp, min = nomin;
    std::atomic<double> minz(std::numeric_limits<double>::max());
    Lpoint min = {0,0.0,0.0,0.0};

    if(isLeaf(qtree))
    {

      if(boxInside2D(*boxMin, *boxMax, qtree)){
        for(Lpoint* p : qtree->points) {
          if (p->z < min.z) {
              min = *p;
          }
          (*numInside)++;
        }
      } else {
        for(Lpoint* p : qtree->points) {
          if(insideBox2D(p, *boxMin, *boxMax))
          {
            if (p->z < min.z) {
                min = *p;
            }
            (*numInside)++;
          }
        }
      }

    } else {

      tbb::parallel_for( 0, static_cast<int>(4) , 1, [&]( int i ) {
        // Check
        bool update_point = true;
        if(boxOverlap2D(*boxMin, *boxMax, qtree->quadrants[i])) {
            Lpoint tmp;
            if(level < 2) // Indicates the parallelization depth
              tmp = findValidMin(qtree->quadrants[i], boxMin, boxMax, numInside, level+1);
            else
              tmp = findValidMin(qtree->quadrants[i], boxMin, boxMax, numInside);

            double local_min = minz.load();
            if (tmp.z < local_min) {
              while(!std::atomic_compare_exchange_weak(&minz, &local_min, tmp.z)){ // Loop until you get it changed
                if(minz.load() < local_min){ // If it fails because someone has a smaller one...
                  update_point = false;
                  break; // Abort!
                }
              }
              if(update_point) min = tmp;
            }
        }

        });

    }

    return min;
}

Lpoint searchNeighborsMinp(Vector2D* point, Qtree qtree, float radius, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMin(qtree, &boxMin, &boxMax, numNeighs, 0);
}

Lpoint searchOverlap(Vector2D* point, Qtree qtree, double radiusX, double radiusY, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radiusX, radiusY, &boxMin, &boxMax);

    return findValidMin(qtree, &boxMin, &boxMax, numNeighs);
}

void countNeighbors(Qtree qtree, Vector2D* boxMin, Vector2D* boxMax, int* numInside)
{
    int i;

    if(isLeaf(qtree))
    {
      if(boxInside2D(*boxMin, *boxMax, qtree)){

          (*numInside) += qtree->points.size();

      } else {
        for(Lpoint* p : qtree->points) {
          if(insideBox2D(p, *boxMin, *boxMax))
          {
            (*numInside)++;
          }
        }
      }

    } else {

        for(i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, qtree->quadrants[i]))
              continue;
            else {
              countNeighbors(qtree->quadrants[i], boxMin, boxMax, numInside);
            }

        }

    }

    return;
}

void countNeighbors2D(Vector2D* point, Qtree qtree, float radius, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    countNeighbors(qtree, &boxMin, &boxMax, numNeighs);

    return;
}



void createQuadrantsF(Qtree qt)
{
    int i = 0;
    Vector2D newCenter;
    float newRadius = qt->radius * 0.5;

    for(i = 0; i < 4; i++)
    {
        newCenter = qt->center;
        newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);

        qt->quadrants[i] = new Qtree_t(newCenter, newRadius);
    }
}


// Insert a point in the qtree creating the appropiate childs
void insertPointF(Lpoint *point, Qtree qtree, float minRadius)
{
    int idx = 0;

    if(isLeaf(qtree))
    {
        if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createQuadrantsF(qtree);
          idx = quadrantIdx(point, qtree);
          insertPointF(point, qtree->quadrants[idx], minRadius);

        } else {
          qtree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      idx = quadrantIdx(point, qtree);
      insertPointF(point, qtree->quadrants[idx], minRadius);
    }
}

void fillQuadrants(Qtree qtree, float minRadius)
{
    int idx;

    for(Lpoint* p : qtree->points)
    {
      idx = quadrantIdx(p, qtree);
      insertPointF(p, qtree->quadrants[idx], minRadius);
    }
    qtree->points.clear();
}

// Insert a point in the qtree creating the appropiate childs
void insertPointF2(Lpoint *point, Qtree qtree, float minRadius, int medSize)
{
    int idx = 0;

    if(isLeaf(qtree))
    {
      if(qtree->points.size() < medSize)
      {
        qtree->points.push_back(point);
      }
      else
      {
        if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createQuadrantsF(qtree);
          fillQuadrants(qtree, minRadius);
          idx = quadrantIdx(point, qtree);
          insertPointF(point, qtree->quadrants[idx], minRadius);

        } else {
          qtree->points.push_back(point);
        }
      }
    }
    else                                // No leaf -> search the correct one
    {
      idx = quadrantIdx(point, qtree);
      insertPointF(point, qtree->quadrants[idx], minRadius);
    }
}

void deleteQtree(Qtree qtree)
{
    int i;
    if(isLeaf(qtree))
    {
      qtree->points.clear();
    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(i = 0; i < 4; i++) {
            // Check
            deleteQtree(qtree->quadrants[i]);
            delete(qtree->quadrants[i]);
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }

    }

    return;
}


void stage1rem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, Qtree qtreeIn, Vector2D min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    Vector2D cellCenter;

    int cellPoints;

    Lpoint newmin = {0,0.0,0.0,0.0};

    Vector2D boxMax, boxMin;

    for(int jj = 0 ; jj < Ccol ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < Crow ; ii++ ){

            cellCenter.x = initX + ii*Displace;

            makeBox(&cellCenter, Wsize*0.5, &boxMin, &boxMax);

            if(insideBox2D(&newmin,boxMin,boxMax)){

              Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

              int old_cellPoints = cellPoints;

              Lpoint tmp = searchOverlap(&oCell, qtreeIn, Displace*0.5, Wsize*0.5, &cellPoints);

              // We're assuming the points were equidistant throughout the cell, which isn't always true.
              cellPoints += (int)(old_cellPoints * Overlap);

              if(tmp.z < newmin.z){
                newmin = tmp;
              }

              // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

            } else {

              newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

            }

            if(cellPoints >= minNumPoints ){
                    minIDs.push_back(newmin.id);
            }
        }

    }

    return;
}

void stage1remCpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, Qtree qtreeIn, Vector2D min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    Vector2D cellCenter;

    int cellPoints;

    Lpoint newmin;

    Vector2D boxMax, boxMin;

    #pragma omp parallel for private(cellCenter,cellPoints,newmin,boxMax,boxMin) schedule(dynamic)
    for(int jj = 0 ; jj < Ccol ; jj++ ){

        cellCenter.x = initX;
        cellCenter.y = initY + jj*Displace;

        newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

        if(cellPoints >= minNumPoints ){
          #pragma omp critical
          {
            minIDs.push_back(newmin.id);
          }
        }

        for(int ii = 1 ; ii < Crow ; ii++ ){

            cellCenter.x = initX + ii*Displace;

            makeBox(&cellCenter, Wsize*0.5, &boxMin, &boxMax);

            if(insideBox2D(&newmin,boxMin,boxMax)){

              Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

              int old_cellPoints = cellPoints;

              Lpoint tmp = searchOverlap(&oCell, qtreeIn, Displace*0.5, Wsize*0.5, &cellPoints);

              // We're assuming the points were equidistant throughout the cell, which isn't always true.
              cellPoints += (int)(old_cellPoints * Overlap);

              if(tmp.z < newmin.z){
                newmin = tmp;
              }

              // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

            } else {

              newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

            }

            if(cellPoints >= minNumPoints ){
              #pragma omp critical
              {
                minIDs.push_back(newmin.id);
              }
            }
        }

    }

    return;
}

void stage1cpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, Qtree qtreeIn, Vector2D min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    Vector2D cellCenter;

    Lpoint newmin;

    int cellPoints;

    #pragma omp parallel for private(cellCenter,cellPoints,newmin) schedule(dynamic)
    for(int jj = 0 ; jj < Ccol ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < Crow ; ii++ ){

            cellCenter.x = initX + ii*Displace;

            newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

            if(cellPoints >= minNumPoints ){

                #pragma omp critical
                {
                    minIDs.push_back(newmin.id);
                }

            }
        }

    }

    return;
}


std::vector<int> stage1tbb(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
    tbb::concurrent_vector<int> v;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // tbb::blocked_range<int>::size_type gs = 4;

    // tbb::affinity_partitioner aff_p;

    tbb::parallel_for( tbb::blocked_range2d<int,int>{0,Ccol,0,Crow},
                       [&](tbb::blocked_range2d<int,int> r ) {

        Lpoint newmin;
        Vector2D cellCenter;
        int cellPoints;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                cellCenter.x = initX + ii*Displace;

                newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

                if( cellPoints >= minNumPoints ){

                    v.push_back(newmin.id);

                }
            }
        }

    });

    return std::vector<int>(v.begin(), v.end());
}

std::vector<int> stage1tbbRem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
    tbb::concurrent_vector<int> v;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // tbb::blocked_range<int>::size_type gs = 4;

    // tbb::affinity_partitioner aff_p;

    tbb::parallel_for( tbb::blocked_range2d<int,int>{0,Ccol,0,Crow},
                       [&](tbb::blocked_range2d<int,int> r ) {

        Lpoint newmin = {0,0.0,0.0,0.0};
        Vector2D cellCenter;
        Vector2D boxMax, boxMin;
        int cellPoints = 0;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

              cellCenter.x = initX + ii*Displace;

              makeBox(&cellCenter, Wsize*0.5, &boxMin, &boxMax);

              if(insideBox2D(&newmin,boxMin,boxMax)){

                Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

                int old_cellPoints = cellPoints;

                Lpoint tmp = searchOverlap(&oCell, qtreeIn, Displace*0.5, Wsize*0.5, &cellPoints);

                // We're assuming the points were equidistant throughout the cell, which isn't always true.
                cellPoints += (int)(old_cellPoints * Overlap);

                if(tmp.z < newmin.z){
                  newmin = tmp;
                }

                // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

              } else {

                newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

              }

              if( cellPoints >= minNumPoints ){

                  v.push_back(newmin.id);

              }

            }
        }

    });

    return std::vector<int>(v.begin(), v.end());
}

std::vector<int> stage1tbb2(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
    tbb::concurrent_vector<int> v;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    std::vector<int> values(Ccol);
    std::iota(std::begin(values), std::end(values), 0);

    tbb::parallel_for_each(values, [&](int jj) {

        // tbb::task_arena nested;

        // tbb::this_task_arena::isolate([&]{

        // nested.execute( [&]{

            for(int ii = 0 ; ii < Crow ; ii++ ){

                Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

                int cellPoints;

                Lpoint newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

                if( cellPoints >= minNumPoints ){

                    v.push_back(newmin.id);

                }
            }
        // });
    });

    return std::vector<int>(v.begin(), v.end());
}

std::vector<int> stage1tg(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
    tbb::concurrent_vector<int> v;
    // std::vector<int> v;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    tbb::task_group tg;

    for(int jj=0; jj<Ccol; jj++){

        tg.run([&,jj]() {

            for(int ii = 0 ; ii < Crow ; ii++ ){

                Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

                int cellPoints;

                Lpoint newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
                // Lpoint newmin = searchNeighborsMinp(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

                if( cellPoints >= minNumPoints ){
                    // std::lock_guard<std::mutex> guard(pool_mutex);
                    v.push_back(newmin.id);

                }
            }

        });
    }

    tg.wait();

    return std::vector<int>(v.begin(), v.end());
}

std::vector<int> stage1tgRem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
    tbb::concurrent_vector<int> v;
    // std::vector<int> v;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    tbb::task_group tg;

    for(int jj=0; jj<Ccol; jj++){

        tg.run([&,jj]() {

            Lpoint newmin = {0,0.0,0.0,0.0};
            Vector2D boxMax, boxMin;
            int cellPoints = 0;

            for(int ii = 0 ; ii < Crow ; ii++ ){

                Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

                makeBox(&cellCenter, Wsize*0.5, &boxMin, &boxMax);

                if(insideBox2D(&newmin,boxMin,boxMax)){

                  Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

                  int old_cellPoints = cellPoints;

                  Lpoint tmp = searchOverlap(&oCell, qtreeIn, Displace*0.5, Wsize*0.5, &cellPoints);

                  // We're assuming the points were equidistant throughout the cell, which isn't always true.
                  cellPoints += (int)(old_cellPoints * Overlap);

                  if(tmp.z < newmin.z){
                    newmin = tmp;
                  }

                  // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

                } else {

                  newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

                }

                if( cellPoints >= minNumPoints ){
                    // std::lock_guard<std::mutex> guard(pool_mutex);
                    v.push_back(newmin.id);

                }
            }

        });
    }

    tg.wait();

    return std::vector<int>(v.begin(), v.end());
}

std::vector<int> stage1tg2(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
    tbb::concurrent_vector<int> v;
    // std::vector<int> v;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    tbb::task_group tg;

    for(int jj=0; jj<Ccol; jj++){

        for(int ii = 0 ; ii < Crow ; ii++ ){

          tg.run([&,ii,jj]() {

            Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

            int cellPoints;

            Lpoint newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

            if( cellPoints >= minNumPoints ){
                // std::lock_guard<std::mutex> guard(pool_mutex);
                v.push_back(newmin.id);

            }

          });

        }

    }

    tg.wait();

    return std::vector<int>(v.begin(), v.end());
}


// std::vector<int> stage1task(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
//   unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
// {
//     tbb::concurrent_vector<int> v;
//     // std::vector<tbb::internal::atomic<int>> counters(Ccol*Crow, 2);
//     //
//     // for (int i=0; i<Crow; i++){
//     //   counters[i]=1;
//     // }
//     // for (int j=0; j<Ccol; j++){
//     //   counters[j*Crow]=1;
//     // }
//     // counters[0] = 0;

//     // for(int j=0; j<Ccol; j++)
//     //   for (int i=0; i<Crow; i++)
//     //     if (i == 0 || j==0) {
//     //       counters[j*Crow+i]=1;
//     //     }
//     // counters[0] = 0;

//     // double Displace = round2d(Wsize*(1-Overlap));

//     double initX = min.x - Wsize/2 + Displace;
//     double initY = min.y - Wsize/2 + Displace;

//     // Vector2D cellCenter = {initX , initY};
//     //
//     // tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root())
//     //                            Cell{0,0,cellCenter,counters,v});

//     // tbb::parallel_for(0, static_cast<int>(Ccol), 1, [&](int j){
//     //   Vector2D cellCenter = {initX , initY + j*Displace};
//     //   tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root())
//     //                              Cell{0,minNumPoints,Ccol,Crow,Wsize,Displace,cellCenter,qtreeIn,v});
//     // });

//     // task_t origin(0,0);
//   	// tbb::parallel_do(&origin, &origin+1, Cell(Ccol,Crow,Wsize,Displace,cellCenter,qtreeIn,v));

//     tbb::task_group tg;

//     for(int jj=0; jj<Ccol; jj++){

//         tg.run([&,jj]() {

//           Vector2D cellCenter = {initX , initY + jj*Displace};

//           tbb::task::spawn_root_and_wait(*new(tbb::task::allocate_root())
//             // Cell{0,minNumPoints,Ccol,Crow,Wsize,Displace,cellCenter,qtreeIn,v});
//             Cell{0,cellCenter,v});

//         });
//     }

//     tg.wait();

//     return std::vector<int>(v.begin(), v.end());
// }



std::vector<int> stage1reduce(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
                        unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  std::vector<int> init;

  return tbb::parallel_reduce(
      tbb::blocked_range2d<int,int>{0,Ccol,0,Crow},
      init,
      [&](const tbb::blocked_range2d<int,int>& r, std::vector<int> v) {

        Lpoint newmin;
        Vector2D cellCenter;
        int cellPoints;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                cellCenter.x = initX + ii*Displace;

                newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

                if( cellPoints >= minNumPoints ){
                    v.push_back(newmin.id);
                }
            }
        }
        return v;
      },
      [](std::vector<int> a, std::vector<int> b) {
         a.insert(a.end(), b.begin(), b.end());
         return a;
      }
  );

}

std::vector<int> stage1reduceRem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
                        unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  std::vector<int> init;

  return tbb::parallel_reduce(
      tbb::blocked_range2d<int,int>{0,Ccol,0,Crow},
      init,
      [&](const tbb::blocked_range2d<int,int>& r, std::vector<int> v) {

        Lpoint newmin = {0,0.0,0.0,0.0};
        Vector2D cellCenter;
        Vector2D boxMax, boxMin;
        int cellPoints = 0;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                cellCenter.x = initX + ii*Displace;

                makeBox(&cellCenter, Wsize*0.5, &boxMin, &boxMax);

                if(insideBox2D(&newmin,boxMin,boxMax)){

                  Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

                  int old_cellPoints = cellPoints;

                  Lpoint tmp = searchOverlap(&oCell, qtreeIn, Displace*0.5, Wsize*0.5, &cellPoints);

                  // We're assuming the points were equidistant throughout the cell, which isn't always true.
                  cellPoints += (int)(old_cellPoints * Overlap);

                  if(tmp.z < newmin.z){
                    newmin = tmp;
                  }

                } else {

                  newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

                }

                if( cellPoints >= minNumPoints ){
                    v.push_back(newmin.id);
                }
            }
        }
        return v;
      },
      [](std::vector<int> a, std::vector<int> b) {
         a.insert(a.end(), b.begin(), b.end());
         return a;
      }
  );

}




unsigned int stage2(unsigned int countMin, std::vector<int>& minIDs){

    unsigned int index = 0;

    int ii,jj,id;

    for( ii=0 ; ii<countMin ; ii=jj ){
        id = minIDs[ii];

        for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

        /* esta es la forma de descartar las posiciones no utilizadas, inicializadas a -1 */
        if(jj-ii > 1){
            minIDs[index]=id;
            index++;
        }
    }

    return index;
}

// std::vector<int> stage2cpp(unsigned int countMin, std::vector<int>& minIDs){

//     unsigned int index = 0;

//     std::vector<int> v;

//     int ii,jj,id;

//     for( ii=0 ; ii<countMin ; ii=jj ){
//         id = minIDs[ii];

//         for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

//         if(jj-ii > 1){
//             v.push_back(id);
//         }
//     }

//     return v;
// }


void stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          std::vector<int>& minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min){

    int cellPoints;

    Vector2D cellCenter;

    Lpoint newmin;

    double initX = min.x + Bsize/2;
    double initY = min.y + Bsize/2;

    for(int jj = 0 ; jj < Ccol ; jj++ ){

      cellCenter.y = initY + jj*Bsize;

      for(int ii = 0 ; ii < Crow ; ii++ ){

        cellCenter.x = initX + ii*Bsize;

        countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);

        // Tengo que hacerlo porque el algoritmo no lo hace
        if(cellPoints == 0){

          newmin = searchNeighborsMin(&cellCenter, qtreeIn, Bsize/2, &cellPoints);

          if(cellPoints>0){
            minGridIDs.push_back(newmin.id);
          }
        }
      }
    }

    return;
}



void stage3cpp(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          std::vector<int>& minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min){

    int cellPoints;

    Vector2D cellCenter;

    Lpoint newmin;

    double initX = min.x + Bsize/2;
    double initY = min.y + Bsize/2;

    #pragma omp parallel for private(cellCenter,cellPoints,newmin) schedule(dynamic)
    for(int jj = 0 ; jj < Ccol ; jj++ ){

      cellCenter.y = initY + jj*Bsize;

      for(int ii = 0 ; ii < Crow ; ii++ ){

          cellCenter.x = initX + ii*Bsize;

          countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);

          // Tengo que hacerlo porque el algoritmo no lo hace
          if(cellPoints == 0){

              newmin = searchNeighborsMin(&cellCenter, qtreeIn, Bsize/2, &cellPoints);

              if(cellPoints>0){

                  #pragma omp critical
                  {
                      minGridIDs.push_back(newmin.id);
                  }

              }
          }
      }
    }

    return;
}

// unsigned int stage3tbb(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
//           std::vector<int>& minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min){

//     // unsigned int addMin = 0;

//     // int cellPoints;

//     // Lpoint cellCenter = {0,0.0,0.0,0.0};

//     // Lpoint** neighbors = NULL;

//     // Lpoint** minimos = NULL;

//     // Lpoint newmin;

//     // int ii,jj;

//     std::atomic<uint32_t> newsum(0);

//     // Qtree mypointer;

//     // #pragma omp parallel for private(cellCenter,cellPoints,newmin) schedule(static)
//     // #pragma omp parallel for private(cellCenter,cellPoints,newmin,mypointer) schedule(dynamic) reduction(+:addMin)
//     // for(int jj = 0 ; jj < Ccol ; jj++ ){
//     tbb::parallel_for(0, static_cast<int>(Ccol),
//         [&](int jj) {

//       Lpoint newmin;
//       Vector2D cellCenter;
//       int cellPoints;
//       // mypointer = qtreeIn;

//       cellCenter.y = min.y + Bsize/2 + jj*Bsize;

//       for(int ii = 0 ; ii < Crow ; ii++ ){

//           cellCenter.x = min.x + Bsize/2 + ii*Bsize;

//           // Candidata a VOID porque no necesito el mínimo, solo el número de puntos encontrados.
//           // newmin = searchNeighborsMin(&cellCenter, grid, Bsize/2, &cellPoints);
//           countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
//           // printf("Numero de elementos de la celda: %d\n", cellPoints );
//           // Tengo que hacerlo porque el algoritmo no lo hace
//           if(cellPoints == 0){
//               // printf("    Tengo una celda vacia en la malla\n" );
//               newmin = searchNeighborsMin(&cellCenter, qtreeIn, Bsize/2, &cellPoints);
//               // newmin = searchNeighborsMinParent(&cellCenter, &mypointer, Bsize/2, &cellPoints);
//               if(cellPoints>0){

//                   // idmin = findMin(neighbors, cellPoints);
//                   // #pragma omp critical
//                   // {
//                   //     minGridIDs[addMin] = newmin.id;
//                   //     addMin++;
//                   // }
//                   minGridIDs[jj*Crow+ii] = newmin.id; //los voy guardando por posición de celda; esto hace que tenga que modificar stage2
//                   // addMin++;
//                   newsum++;
//                   // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
//               }
//           }
//       }
//     });

//     return newsum.load();

// }

std::vector<int> stage3reduce(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
                            Qtree qtreeIn, Qtree grid, Vector2D min)
{
  double initX =  min.x + Bsize/2;
  double initY = min.y + Bsize/2;

  std::vector<int> init;

  return tbb::parallel_reduce(
      tbb::blocked_range2d<int,int>{0,Ccol,0,Crow},
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

                countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);

                if( cellPoints == 0 ){

                    newmin = searchNeighborsMin(&cellCenter, qtreeIn, Bsize/2, &cellPoints);

                    if( cellPoints > 0 ){
                        v.push_back(newmin.id);
                    }
                }
            }
        }
        return v;
      },
      [](std::vector<int> a, std::vector<int> b) {
         a.insert(a.end(), b.begin(), b.end());
         return a;
      }
  );

}

int save_file(std::string file_name, std::vector<int>& ids, Lpoint** pointer)
{
  // unsigned int size = pointList.size();
  std::ofstream out;
  out.precision(std::numeric_limits<double>::digits10);
  out.open(file_name);
  // Point aux;
  if(!out.is_open())
    return -1;
  for( int i : ids ){
    out << (*pointer)[i].x << " " << (*pointer)[i].y << " " << (*pointer)[i].z << std::endl;
    // fprintf(out, "%.2f %.2f %.2f\n", (*pointer)[i].x, (*pointer)[i].y, (*pointer)[i].z);
  }
  out.close();
  if(out.is_open())
    return -1;
  return 0;
}

int check_results(std::string filename, std::vector<int>& ids, Lpoint** pointer, float Displace)
{
  std::string line;
  // char delim;
  Vector3D aux;
  unsigned int npoints = 0;
  std::vector<Vector3D> vgold, v;

  std::ifstream input(filename);
  if (input.fail()) {
      std::cout << "The GOLD point list doesn't exits" << std::endl;
      return 0;
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
    return 0;

  std::sort(vgold.begin(), vgold.end(), [](const Vector3D a, const Vector3D b){
    return a.z < b.z;
  });

  // for(auto p1 : v){
  //   printf("%lf %lf %lf\n", p1.x, p1.y, p1.z);
  // }

  for( int i : ids ){
    aux.x = (*pointer)[i].x;
    aux.y = (*pointer)[i].y;
    aux.z = (*pointer)[i].z;
    v.push_back(aux);
  }

  std::sort(v.begin(), v.end(), [](const Vector3D a, const Vector3D b){
    return a.z < b.z;
  });

  int count = 0;
  // float Displace = 1.0;
  Vector2D boxMin, boxMax;
  for(auto p : v){
    Vector2D aux = {p.x, p.y};
    makeBox(&aux, Displace*0.5, &boxMin, &boxMax);
    for(auto pg : vgold){
      if(insideBox2D(&pg,boxMin,boxMax)){
        // if(fabs(p.z - pg.z) < 0.01){
          count++;
          break;
        }
    }
  }

  double rate = count/((double)(npoints))*100.0;
  printf("%d points correct; %.2f%%\n", count, rate);

  return 0;

  // return std::equal(v.begin(), v.end(), v2.begin(), [](const Vector3D a, const Vector3D b){
  //   return fabs(a.z - b.z) < 0.01;
  // });
}

uint64_t read_points(std::string filename, Lpoint** point_cloud) 
{
  
  std::ifstream input(filename);
  if (input.fail())
    return -1;

  std::string line;
  uint64_t id = 0;
  
  while(getline(input,line)) {
    std::stringstream s(line); 
    (*point_cloud)[id].id = id;
    s >> (*point_cloud)[id].x >> (*point_cloud)[id].y >> (*point_cloud)[id].z;
    id++;
  }

  input.close();
  if(input.is_open())
    return -1;

  return id;
}

int readXYZfile(std::string filename, Lpoint* & point_cloud, unsigned int & Npoints, Vector2D &min, Vector2D &max) 
{
  FILE* fileXYZ;
  if((fileXYZ = fopen(filename.c_str(),"r")) == NULL){
    printf("Unable to open file!\n");
    return -1;
  }
  if ( filename.find("Arzua.xyz") != std::string::npos || filename.find("Alcoy.xyz") != std::string::npos || 
       filename.find("BrionF.xyz") != std::string::npos || filename.find("BrionU.xyz") != std::string::npos ){
    printf("Read header...\n");
    if(fscanf(fileXYZ, "%d\n%lf\n%lf\n%lf\n%lf\n",&Npoints, &min.x, &max.x, &min.y, &max.y) < 5){
        printf("Imposible to read header values\n");
        return -1;
    }
  }
  // Allocate memory for the LiDAR points
  point_cloud = (Lpoint*)malloc(Npoints*sizeof(Lpoint));
  if(point_cloud == NULL){
    printf("Error allocating memory for the points\n");
    return -1;
  }
  printf("Reading points...\n");

  for(int i=0; i<Npoints ; i++){
    point_cloud[i].id = i;
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


// void prueba1( uQtree2& qtree )
// {
//   printf("myfloat %g\n", qtree->radius);
//   qtree.reset( new nHoja2({0.0,0.0,0.0}, 33.0) );
//   // printf("myfloat %g\n", (qtree)->radius);
//   // qtree.reset();
//   return;
// }
