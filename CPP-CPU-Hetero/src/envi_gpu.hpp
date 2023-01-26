#ifndef ENVI_GPU_HPP

#define ENVI_GPU_HPP

// #define TBB_PREVIEW_MEMORY_POOL 1
#define TBB_PREVIEW_GLOBAL_CONTROL 1

#include <sycl/sycl.hpp>
#include <iostream>
// #include <limits>
#include <cmath>
// #include <omp.h>
#include <vector>
#include <unistd.h>
// #include <functional>
// #include <algorithm>
#include <string.h>
// #include <numeric>
#include <atomic>
#include <fstream>
#include <sstream>
#include <chrono>
#include <numeric>

// #include <CL/sycl/backend/cuda.hpp>

#include "tbb/global_control.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for.h"
#include "tbb/concurrent_vector.h"
#include "tbb/blocked_range2d.h"
#include "tbb/task_group.h"
// #include "tbb/cache_aligned_allocator.h"
// #include "tbb/scalable_allocator.h"
#include "tbb/tbbmalloc_proxy.h"
// #include "tbb/task_arena.h"
// #include "tbb/memory_pool.h"
// #include "tbb/atomic.h"

// Processes input arguments
#include "../include/cxxopts.hpp"

using namespace sycl;





// NEW TYPES

typedef struct
{
    double x;
    double y;

} Vector2D;


typedef struct
{
    double x;
    double y;
    double z;

} Vector3D;


typedef struct
{
    uint32_t id;
    double x;
    double y;
    double z;

} Lpoint;


template<typename T>
struct addr_t
{

  T first;
  T second;
  T third;
  T fourth;

  T& operator[] (int pos) {
  
  if(pos==0)
    return first;
  else if(pos==1)
    return second;
  else if(pos==2)
    return third;
  else
    return fourth;

  }

};

const usm::alloc global_alloc = usm::alloc::host;
// const usm::alloc global_alloc = usm::alloc::shared;
// const usm::alloc global_alloc = usm::alloc::device;

struct Qtree_t; // CPU node

typedef struct Qtree_t* Qtree;

struct Qtree_t {
  Qtree quadrants[4];
  // std::vector<Qtree> quadrants;
  Qtree parent;
  // Vector3D center;
  Vector2D center;
  float radius;
  // std::vector<Lpoint*> points;
  tbb::concurrent_vector<Lpoint*> points;

  Qtree_t(){};
  Qtree_t( Qtree p, Vector2D& c, float r ) : 
              parent(p), center(c), radius(r) {

      // quadrants.reserve(4);

      for(int i=0; i<4; ++i)
        quadrants[i] = NULL;

    };
};

struct QtreeG4_t; // GPU node

typedef struct QtreeG4_t* QtreeG4;

struct QtreeG4_t {
  // addr_t<QtreeG4> quadrants;
  // std::vector<QtreeG4, usm_allocator<QtreeG4, global_alloc>> quadrants;
  QtreeG4 parent;
  // QtreeG4* quadrants;
  QtreeG4 quadrants;
  Vector2D center;
  float radius;
  // point_vector_t points;
  int numPts;
  // Lpoint** points;
  Lpoint* points;
  Lpoint* min;
  
  QtreeG4_t( QtreeG4 p, Vector2D& c, float r ) : 
              parent(p), center(c), radius(r) {

      quadrants = NULL;
      numPts = 0;
      points = NULL;
      min = NULL;

    };
};

struct QtreeG5_t; // GPU node with indexes

typedef struct QtreeG5_t* QtreeG5;

struct QtreeG5_t {
  // addr_t<QtreeG4> quadrants;
  // std::vector<QtreeG4, usm_allocator<QtreeG4, global_alloc>> quadrants;
  int32_t parent;
  // QtreeG4* quadrants;
  int32_t quadrants;
  Vector2D center;
  float radius;
  // point_vector_t points;
  int numPts;
  // Lpoint** points;
  int32_t points;
  int32_t min;
  
  QtreeG5_t( int32_t p, Vector2D& c, float r ) : 
              parent(p), center(c), radius(r) {

      quadrants = -1;
      numPts = 0;
      points = -1;
      min = -1;

      // for(int i=0; i<4; ++i)
      //   quadrants[i] = NULL;

    };
};





// GLOBAL DEFINITIONS

uint32_t Wsize;
uint32_t Bsize;
uint32_t minNumPoints;
uint32_t nCols;
double initX;
double initY;
double Overlap;
double Displace;
int* minIDs = NULL;

uint16_t is_gpu_used = 0;

uint64_t Npoints = 0;

Lpoint* point_cloud = NULL;

Qtree cpu_qtree = NULL;

// uint64_t cpu_tree_nodes = 0;
std::atomic<uint64_t> cpu_tree_nodes = {0};
// tbb::internal::atomic<uint64_t> cpu_tree_nodes = {0};

// Lpoint** array_all_points = NULL;
Lpoint* array_all_points = NULL;
// static constexpr Lpoint* array_all_points = NULL;

/* qtree for GPU with aligned nodes and points */
QtreeG4 aligned_qtree = NULL;

QtreeG4 array_pointers = NULL;

uint64_t copied_nodes = 0;

uint64_t num_pts_copiados = 0;

/* qtree for GPU with aligned nodes and points, which stores only
  the indexes */
QtreeG5 array_indexes = NULL;

int32_t cuda_copied_nodes = 0;

int32_t cuda_num_pts_copiados = 0;


// gpu_selector selector;
// // cpu_selector hselector;
// //default_selector selector;
// //host_selector selector;

// queue q(selector);
sycl::queue gpu_queue(gpu_selector_v /*, propList*/);

// sycl::device qdevice = gpu_queue.get_device();
// sycl::context qcontext = gpu_queue.get_context();

tbb::task_group global_tg;





// FUNCTIONS

void* mallocWrap(size_t size);

Lpoint gpuSearchNeighborsMin(Vector2D& center, QtreeG4 qtree, float radiusX, float radiusY, int& numNeighs);

template<typename accNode_t, typename accPoint_t>
Lpoint gpuSearchIndex( Vector2D& center, accNode_t aNodes, accPoint_t aPoints, float radiusX, float radiusY, int& numNeighs);


int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int exponent(double value, int base) 
{
  int i = 0;
  while(pow(base,i) < value){
    i++;
  }
  return i;
}

template<typename it_t>
it_t siftIterator(it_t myiterator, int stride) 
{
  it_t newIt = myiterator;
  std::advance(newIt, stride);
  return newIt;
}

double round2d(double z){
  return round(z*100.0)/100.0;
}


/* Makes a box of equal sides with two points */
void makeBox(Vector2D& center, float radius, Vector2D& min, Vector2D& max)
{
    min.x = center.x - radius;
    min.y = center.y - radius;
    // min.z = center.z - radius;

    max.x = center.x + radius;
    max.y = center.y + radius;
    // max.z = point.z + radius;
}


/* Makes a box of different sides with two points */
void makeBox(Vector2D& center, double radiusX, double radiusY, Vector2D& min, Vector2D& max)
{

    min.x = center.x - radiusX;
    min.y = center.y - radiusY;

    max.x = center.x + radiusX;
    max.y = center.y + radiusY;

}


/* Checks if a point is contained in the search area */
int insideBox2D(Vector3D& point, Vector2D& min, Vector2D& max)
{
    if(point.x > min.x && point.y > min.y)
    {
        if(point.x < max.x && point.y < max.y)
        {
            return 1;
        }
    }

    return 0;
}


/* Checks if a point is contained in the search area */
int insideBox2D(Lpoint* point, Vector2D& min, Vector2D& max)
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


/* Checks if the node is partially overlapped by the search area */
template<typename pointer_t>
int boxOverlap2D(Vector2D& boxMin, Vector2D& boxMax, pointer_t qt)
{
    if(qt->center.x + qt->radius < boxMin.x ||
       qt->center.y + qt->radius < boxMin.y)
        return 0;

    if(qt->center.x - qt->radius > boxMax.x ||
       qt->center.y - qt->radius > boxMax.y)
        return 0;

    return 1;
}


/* Checks if the search area is completely overlapped by the node */
template<typename pointer_t>
int boxInside2D(Vector2D& boxMin, Vector2D& boxMax, pointer_t qt)
{
    if(qt->center.x + qt->radius > boxMax.x ||
       qt->center.y + qt->radius > boxMax.y)
        return 0;

    if(qt->center.x - qt->radius < boxMin.x ||
       qt->center.y - qt->radius < boxMin.y)
        return 0;

    return 1;
}


/* Checks if the node is completely overlapped by the search area */
template<typename pointer_t>
int boxTotalOverlap2D(Vector2D& boxMin, Vector2D& boxMax, pointer_t qt)
{
    if(qt->center.x + qt->radius < boxMax.x ||
       qt->center.y + qt->radius < boxMax.y)
        return 0;

    if(qt->center.x - qt->radius > boxMin.x ||
       qt->center.y - qt->radius > boxMin.y)
        return 0;

    return 1;
}



template<typename pointer_t>
int isLeaf(pointer_t qt)
{
    return qt->quadrants[0] == NULL;
}


/** Calculate the radius in each axis and save the max radius of the bounding box */
Vector2D getRadius(Vector2D min, Vector2D max, float *maxRadius)
{
    Vector2D radii;

    radii.x = (max.x - min.x) / 2.0;
    radii.y = (max.y - min.y) / 2.0;
    // radii.z = (max.z - min.z) / 2.0;

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
    // center.z = min.z + radius.z;

    return center;
}


/* Identifies the index of the quadrant where the point should
 be inserted */
template<typename pointer_t>
int quadrantIdx(Lpoint *point, pointer_t qtree)
{
    int child = 0;

    if(point->x >= qtree->center.x) child |= 2;
    if(point->y >= qtree->center.y) child |= 1;
    // if(point->z >= qtree->center.z) child |= 1;

    return child;
}


/* Recursive search of the minimum in an area */
template<typename pointer_t>
Lpoint findValidMin(pointer_t qtree, Vector2D& boxMin, Vector2D& boxMax, int& numInside)
{
    // Lpoint tmp, min = nomin;
    Lpoint tmp, min = {0,0.0,0.0,99999.0};

    if(isLeaf(qtree))
    {

      if(boxInside2D(boxMin, boxMax, qtree)){
        for(Lpoint* p : qtree->points) {
          if (p->z < min.z) {
              min = *p;
          }
          numInside++;
        }
      } else {
        for(Lpoint* p : qtree->points) {
          if(insideBox2D(p, boxMin, boxMax))
          {
            if (p->z < min.z) {
                min = *p;
            }
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
                tmp = findValidMin(qtree->quadrants[i], boxMin, boxMax, numInside);
                if (tmp.z < min.z) {
                    min = tmp;
                }
            }

        }

    }

    return min;
}

template<typename pointer_t>
Lpoint searchNeighborsMin(Vector2D& center, pointer_t qtree, float radius, int& numNeighs)
{
    Vector2D boxMin, boxMax;

    numNeighs = 0;
    makeBox(center, radius, boxMin, boxMax);

    return findValidMin(qtree, boxMin, boxMax, numNeighs);
}

/* Same as searchNeighborsMin, but creating a search area with
 uneven sides for the algorithm where we remember the minimum */
Lpoint searchOverlap(Vector2D& center, Qtree qtree, double radiusX, double radiusY, int& numNeighs)
{
    Vector2D boxMin, boxMax;

    numNeighs = 0;
    makeBox(center, radiusX, radiusY, boxMin, boxMax);

    return findValidMin(qtree, boxMin, boxMax, numNeighs);
}


/* Recursive function to count the number of points inside
 a search area*/
void countNeighbors(Qtree qtree, Vector2D& boxMin, Vector2D& boxMax, int& numInside)
{
    int i;

    if(isLeaf(qtree))
    {
      if(boxInside2D(boxMin, boxMax, qtree)){

          numInside += qtree->points.size();

      } else {
        for(Lpoint* p : qtree->points) {
          if(insideBox2D(p, boxMin, boxMax))
          {
            numInside++;
          }
        }
      }

    } else {

        for(i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(boxMin, boxMax, qtree->quadrants[i]))
              continue;
            else {
              countNeighbors(qtree->quadrants[i], boxMin, boxMax, numInside);
            }

        }

    }

    return;
}

void countNeighbors2D(Vector2D& center, Qtree qtree, float radius, int& numNeighs)
{
    Vector2D boxMin, boxMax;

    numNeighs = 0;
    makeBox(center, radius, boxMin, boxMax);

    countNeighbors(qtree, boxMin, boxMax, numNeighs);

    return;
}


/* Recursive function to count the number of points in the leaf nodes */
void register_leafs(Qtree qtree, Vector2D& boxMin, Vector2D& boxMax, std::vector<int>& v)
{
    int i;

    if(isLeaf(qtree))
    {

      v.push_back(qtree->points.size());

    } else {

        for(i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(boxMin, boxMax, qtree->quadrants[i]))
              continue;
            else {
              register_leafs(qtree->quadrants[i], boxMin, boxMax, v);
            }

        }

    }

    return;
}

void register_data(Vector2D& center, Qtree qtree, float radius, std::vector<int>& v)
{
    Vector2D boxMin, boxMax;

    makeBox(center, radius, boxMin, boxMax);

    register_leafs(qtree, boxMin, boxMax, v);

    return;
}


/* Initializes one node */
Qtree createQtree(Qtree parent, Vector2D center, float radius)
{
    Qtree qt = new Qtree_t;
    // Qtree qt = node_allocator.allocate(1);

    qt->center = center;
    qt->radius = radius;
    qt->parent = parent;

    // qt->quadrants.reserve(4);
    
    for( int i = 0; i < 4; i++)
      qt->quadrants[i] = NULL;

    cpu_tree_nodes++;

    return qt;
}


/* Initializes the four quadrants */
void createQuadrants(Qtree qt)
{
    Vector2D newCenter;
    float newRadius = qt->radius * 0.5;

    for( int i = 0; i < 4; i++)
    {
        newCenter = qt->center;
        newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);

        // qt->quadrants[i] = createQtree(qt, newCenter, newRadius);
        qt->quadrants[i] = new Qtree_t(qt, newCenter, newRadius);

        cpu_tree_nodes++;

    }
}



void insertPoint(Lpoint* point, Qtree qtree, float minRadius);

/* Re-inserts the points when its exceed the set max. size  */
void fillQuadrants(Qtree qtree, float minRadius)
{

    for(Lpoint* p : qtree->points)
    {
      int idx = quadrantIdx(p, qtree);
      insertPoint(p, qtree->quadrants[idx], minRadius);
    }

    qtree->points.clear();
}

/* Inserts a point in the qtree creating the appropiate childs
 by setting a minimum area*/
void insertPoint(Lpoint* point, Qtree qtree, float minRadius)
{
    int idx = 0;

    if(isLeaf(qtree))
    {
        // printf("octante hoja nivel %d\n",nivel);
        if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createQuadrants(qtree);
          // fillOctants(qtree);
          idx = quadrantIdx(point, qtree);
          insertPoint(point, qtree->quadrants[idx], minRadius);

        } else {
          qtree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
      idx = quadrantIdx(point, qtree);
      insertPoint(point, qtree->quadrants[idx], minRadius);
    }
}


void insertPoint2(Lpoint *point, Qtree qtree, int maxNumber);

/* Re-inserts the points when its exceed the set max. size  */
void fillQuadrants(Qtree qtree, int maxNumber)
{

    for(Lpoint* p : qtree->points)
    {
      int idx = quadrantIdx(p, qtree);
      insertPoint2(p, qtree->quadrants[idx], maxNumber);
    }

    qtree->points.clear();
}


/* Inserts a point in the qtree creating the appropiate childs
 by setting a maximum number of points*/
void insertPoint2(Lpoint* point, Qtree qtree, int maxNumber)
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
        fillQuadrants(qtree, maxNumber);
        int idx = quadrantIdx(point, qtree);
        insertPoint2(point, qtree->quadrants[idx], maxNumber);

      }
    }
    else                                // No leaf -> search the correct one
    {
      int idx = quadrantIdx(point, qtree);
      insertPoint2(point, qtree->quadrants[idx], maxNumber);
    }
}


/* Creates nodes up to "max level" and from this level it starts
saving the pointer to the nodes in a vector that later on will be processed in parallel */
void insert_leafs(Qtree qtree, int maxlevel, int level, std::vector<Qtree>& n_work)
{

    if(level < maxlevel)
    {
      createQuadrants(qtree);

      for(int i = 0; i < 4; i++)
        insert_leafs(qtree->quadrants[i], maxlevel, level+1, n_work);

    }
    else
    { 
      n_work.push_back(qtree);
    }

}


/* Creates the tree in parallel using parallel_for */
Qtree parallel_qtree_creation( int level, Vector2D center, float radius, int maxNumber )
{

  Qtree root = createQtree( NULL, center, radius);

  std::vector<Qtree> n_work;

  insert_leafs( root, level, 0, n_work );

  printf("    N tareas: %zu\n", n_work.size());

  tbb::parallel_for( 0, static_cast<int>(n_work.size()), 1,
    [&](int id){

    Qtree qt = n_work[id];

    Vector2D boxMin, boxMax;
    makeBox(qt->center, qt->radius, boxMin, boxMax);

    for(int i = 0; i < Npoints; i++) {

      Lpoint* point = &point_cloud[i];

      if( insideBox2D(point, boxMin, boxMax) ) {
        // insertPoint(point, qt, 0.5);
        insertPoint2(point, qt, maxNumber);
      }

    }

  });

  return root;
}


/* Accumulates the points on a certain level of the tree */
void insertPoint3(Lpoint* point, Qtree qtree)
{
    if(isLeaf(qtree))
    {
      // std::cout << "  Level: " << level << std::endl << std::flush;
      qtree->points.push_back(point);
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
      int idx = quadrantIdx(point, qtree);
      insertPoint3(point, qtree->quadrants[idx]);
    }
}

//This function can be called with node_delimiter as MinRadius or MaxNumber
template<typename delimiter_t> // int for MaxNumber or float for MinRadius
Qtree parallel_qtree_pf2( int level, Vector2D center, float radius, delimiter_t node_delimiter )
{

  Qtree root = createQtree( NULL, center, radius);

  std::vector<Qtree> n_work;
//create the "transitory" leaves up to tree-level "level" (sequential) and store these leaves in the n_work vector
  insert_leafs( root, level, 0, n_work );

  // std::cout << "  N tareas: " << n_work.size() << std::endl << std::flush;
//traverse the LiDAR points in parallel in the transitory leaves
  tbb:: parallel_for( tbb::blocked_range<int>{0, static_cast<int>(Npoints)},
                      [root](tbb::blocked_range<int> r ) {

    int end = r.end();

    for(int i = static_cast<int>(r.begin()); i < end; i++) {
      insertPoint3(&point_cloud[i], root);
    }

  });

//Finally, traverse in parallel the transitory leaves and finish up the tree
//that hangs from them storing the points in the final leaves
  tbb::parallel_for( 0, static_cast<int>(n_work.size()), 1,
    [&](int id){

    Qtree qt = n_work[id];

    createQuadrants(qt);
    //Depending on the type of node delimiter it will call to the MinRadius or MaxNumber version
    fillQuadrants(qt, node_delimiter);

  });

  return root;
}


template<typename delimiter_t>
Qtree parallel_qtree_stage3( int level, Vector2D center, float radius, delimiter_t node_delimiter, int* minIDs, int countGPU )
{

  Qtree root = createQtree( NULL, center, radius);

  std::vector<Qtree> n_work;

  insert_leafs( root, level, 0, n_work );

  // std::cout << "  N tareas: " << n_work.size() << std::endl << std::flush;

  tbb:: parallel_for( tbb::blocked_range<int>{0, static_cast<int>(countGPU)},
                      [root, minIDs](tbb::blocked_range<int> r ) {

    int end = r.end();

    for(int i = static_cast<int>(r.begin()); i < end; i++) {
      insertPoint3(&point_cloud[minIDs[i]], root);
    }

  });


  tbb::parallel_for( 0, static_cast<int>(n_work.size()), 1,
    [&](int id){

    Qtree qt = n_work[id];

    createQuadrants(qt);
    fillQuadrants(qt, node_delimiter);

  });

  return root;
}


void go_deeper(Qtree qt, int maxNumber, int limit);

Qtree parallel_qtree_pf3( int level, Vector2D center, float radius, int maxNumber, int limit )
{

  Qtree root = createQtree( NULL, center, radius);

  std::vector<Qtree> n_work;

  insert_leafs( root, level, 0, n_work );

  // std::cout << "  N tareas: " << n_work.size() << std::endl << std::flush;

  tbb:: parallel_for( tbb::blocked_range<int>{0, static_cast<int>(Npoints)},
                      [root](tbb::blocked_range<int> r ) {

    int end = r.end();

    for(int i = static_cast<int>(r.begin()); i < end; i++) {
      insertPoint3(&point_cloud[i], root);
    }

  });


  tbb::parallel_for( 0, static_cast<int>(n_work.size()), 1,
    [&n_work,maxNumber,limit](int id){

    Qtree qt = n_work[id];

    createQuadrants(qt);

    if(qt->points.size() < limit) {
      fillQuadrants(qt, maxNumber);
    }
    else {
      go_deeper(qt, maxNumber, limit);
      qt->points.clear();
    }

  });

  return root;
}

void go_deeper(Qtree qt, int maxNumber, int limit)
{
  std::vector<Qtree> n_work;

  insert_leafs( qt, 2, 0, n_work );

  // Lpoint** points = qt->points.data(); 
  // auto& points = qt->points; 

  // tbb:: parallel_for( tbb::blocked_range<int>{0, static_cast<int>(size)},
  //                     [qt](tbb::blocked_range<int> r ) {

  //   int end = r.end();

  //   for(int i = static_cast<int>(r.begin()); i < end; i++) {
  //     insertPoint3(qt->points[i], qt, 0);
  //   }

  // });

  tbb:: parallel_for( 
    tbb::blocked_range< tbb::concurrent_vector<Lpoint*>::iterator >{qt->points.begin(),qt->points.end()},
                      [qt]( const tbb::blocked_range< tbb::concurrent_vector<Lpoint*>::iterator >& r ) {

    for(auto& p : r) {
      insertPoint3(p, qt);
    }

  });

  tbb::parallel_for( 0, static_cast<int>(n_work.size()), 1,
    [&n_work,maxNumber,limit](int id){

    Qtree leaf = n_work[id];

    createQuadrants(leaf);

    if(leaf->points.size() < limit) {
      fillQuadrants(leaf, maxNumber);
    }
    else {
      go_deeper(leaf, maxNumber, limit);
      leaf->points.clear();
    }

  });

  return;
}


/* Creates nodes up to "max level" and form this level it starts
creating tasks to insert points */
void insert_leafs(Qtree qtree, int maxlevel, int level, int maxNumber)
{

    if(level < maxlevel)
    {
      createQuadrants(qtree);

      for(int i = 0; i < 4; i++)
        insert_leafs(qtree->quadrants[i], maxlevel, level+1, maxNumber);

    }
    else
    {

      Lpoint* captured_cloud = point_cloud;
      uint64_t size = Npoints;

      global_tg.run( [qtree,captured_cloud,size,maxNumber]{

        Vector2D boxMin, boxMax;
        makeBox(qtree->center, qtree->radius, boxMin, boxMax);

        for(int i = 0; i < size; i++) {

          Lpoint* p = &captured_cloud[i];

          if( insideBox2D(p, boxMin, boxMax) ) {
            // insertPoint(point, qtree, 0.5);
            insertPoint2(p, qtree, maxNumber);
          }

        } //for

      }); //task
      
    } //else
      
}

/* Creates the tree in parallel using task_group */
Qtree parallel_qtree_creationtg( int level, Vector2D center, float radius, int maxNumber )
{

  Qtree root = createQtree( NULL, center, radius);

  insert_leafs( root, level, 0, maxNumber);

  global_tg.wait();

  return root;
}


/* Creates nodes up to "max level" and form this level it starts
creating tasks that first, look for all the points in that area,
and then start the insertion creating the tree */
void search_and_save(Qtree qtree, int maxlevel, int level, int maxNumber)
{

    if(level < maxlevel)
    {
      createQuadrants(qtree);

      for(int i = 0; i < 4; i++)
        search_and_save(qtree->quadrants[i], maxlevel, level+1, maxNumber);

    }
    else
    {

      Lpoint* captured_cloud = point_cloud;
      uint64_t size = Npoints;

      global_tg.run( [qtree, captured_cloud, size, maxNumber]{

        Vector2D boxMin, boxMax;
        makeBox(qtree->center, qtree->radius, boxMin, boxMax);

        for(int i = 0; i < size; i++) {

          Lpoint* p = &captured_cloud[i];

          if( insideBox2D(p, boxMin, boxMax) ) {
            // insertPoint(point, qtree, 0.5);
            // insertPoint2(p, qtree, maxNumber);
            qtree->points.push_back(p);
          }

        } //for

        createQuadrants(qtree);
        fillQuadrants(qtree, maxNumber);

      }); //task
      
    } //else
      
}

Qtree parallel_qtree_creationtg2( Vector2D center, float radius, int maxlevel, int maxNumber )
{

  Qtree root = new Qtree_t( NULL, center, radius);

  search_and_save( root, maxlevel, 0, maxNumber);

  global_tg.wait();

  return root;
}


/* Recursively deletes the tree created by CPU */
void deleteQtree(Qtree qtree)
{
    if(isLeaf(qtree))
    {

      qtree->points.clear();
      // qtree->quadrants.clear();

    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; i++) {
            // Check
            deleteQtree(qtree->quadrants[i]);
            delete(qtree->quadrants[i]);
        }
        // qtree->quadrants.clear();
    }

    return;
}


/*Este se utiliza cuando NO tengo los cuadrantes consecutivos en memora*/
/* It looks for the posisition of the current node in the 
quadrant array of the parent node */
// void findPosition(QtreeG4 parent, QtreeG4& current, int& idx) 
// {

//   if(parent->quadrants[0] == current){
//     idx=1;
//     current = parent->quadrants[1];
//   }
//   else if(parent->quadrants[1] == current){
//     idx=2;
//     current = parent->quadrants[2];
//   }
//   else if(parent->quadrants[2] == current){
//     idx=3;
//     current = parent->quadrants[3];
//   } 
//   else {
//     idx = 4; // termino, pero dejo listo current
//     current = parent->quadrants[3];
//   }

// }


/* Este se utiliza cuando tengo los cuadrantes consecutivos en memoria */
/* It looks for the posisition of the current node in the 
quadrant array of the parent node */
void findPosition(QtreeG4 parent, QtreeG4& current, int& idx) 
{

  if(parent->quadrants == current){
    // std::advance(iterador, 1);
    idx=1;
    current = &(parent->quadrants[1]);
  }
  else if(&(parent->quadrants[1]) == current){
    // std::advance(iterador, 2);
    idx=2;
    current = &(parent->quadrants[2]);
  }
  else if(&(parent->quadrants[2]) == current){
    // std::advance(iterador, 3);
    idx=3;
    current = &(parent->quadrants[3]);
  } 
  else {
    //lo pongo apuntando al último
    // std::advance(iterador, 3);
    idx = 4; // termino, pero dejo listo current
    current = &(parent->quadrants[3]);
  }

  // int child = 0;

  // Vector2D point = current->center;

  // if(point.x >= parent->center.x) child |= 2;
  // if(point.y >= parent->center.y) child |= 1;

  // idx = child + 1;

  // if(idx < 4)
  //   current = &(parent->quadrants[idx]);
}


/* Sequential Stage1 with mins vector*/
template<typename pointer_t>
void stage1s(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, std::vector<int>& minIDs, pointer_t qtreeIn, Vector2D min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    Vector2D cellCenter;

    Lpoint newmin;

    int cellPoints;

    for(int jj = 0 ; jj < nRows ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < nCols ; ii++ ){

          // printf("CELDA %d,%d <-------------------------------------------------------\n",jj,ii);

            cellCenter.x = initX + ii*Displace;

            // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
            // newmin = cpuSearchNeighborsMin(cellCenter, qtreeIn, Wsize/2, cellPoints);
            newmin = gpuSearchNeighborsMin(cellCenter, qtreeIn, Wsize*0.5, Wsize*0.5, cellPoints);
            // newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

            if(cellPoints >= minNumPoints ){
              // printf("minimo: %d", newmin.id);
              minIDs.push_back(newmin.id);
            }
        }

    }

    return;
}


/* Parallel stage1 with parallel_for that splits the work in blocks of cells (1D),
 instead of in blocks of rows (2D)*/
template<typename pointer_t>
void stage1tbbOne(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, int* minIDs, pointer_t qtreeIn, Vector2D min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    int cellEnd = nRows*nCols;

    // tbb::parallel_for( 0, cellEnd, [&]( int i ) {
    tbb:: parallel_for( tbb::blocked_range<int>{0,cellEnd},
                      [&](tbb::blocked_range<int> r ) {

      int cellPoints;

      int end = static_cast<int>(r.end());
      
      for(int i = static_cast<int>(r.begin()); i<end; i++ ) {

        Vector2D cellCenter = {initX + (i%nCols)*Displace, initY + (int)(i/nCols)*Displace};

        // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
        // newmin = cpuSearchNeighborsMin(cellCenter, qtreeIn, Wsize/2, cellPoints);
        Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtreeIn, Wsize*0.5, Wsize*0.5, cellPoints);
        // newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

        if(cellPoints >= minNumPoints ){
          // printf("minimo: %d", newmin.id);
          minIDs[i] = newmin.id;
        }

      }

    });

    return;
}


/* Heterogeneous stage1 that splits the work in blocks of cells (1D),
 instead of in blocks of rows (2D)*/
template<typename pointer_t>
void stage1heterOne(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, int* minIDs, pointer_t qtree, Vector2D min, float rate)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    is_gpu_used = 1;

    int cellEnd = nRows*nCols;

    tbb::task_group tg;

    int chunk = static_cast<int>(cellEnd * rate);

    tg.run( [&]{

      gpu_queue.parallel_for(range<1>(chunk), [=](sycl::id<1> id) { 

        int i = static_cast<int>(id[0]);

        int cellPoints;

        Vector2D cellCenter = {initX + (i%nCols)*Displace, initY + (int)(i/nCols)*Displace};

        Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);

        if(cellPoints >= minNumPoints ){
          minIDs[i] = newmin.id;
        }

      }).wait();

    });

    tg.run( [&]{

      // tbb::parallel_for( chunk, cellEnd, [&]( int i ) {
      tbb:: parallel_for( tbb::blocked_range<int>{chunk,cellEnd},
                       [&](tbb::blocked_range<int> r ) {

        int end = static_cast<int>(r.end());
        
        for(int i = static_cast<int>(r.begin()); i<end; i++ ) {

          int cellPoints;

          Vector2D cellCenter = {initX + (i%nCols)*Displace, initY + (int)(i/nCols)*Displace};

          Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);

          if(cellPoints >= minNumPoints ){

            minIDs[i] = newmin.id;

          }

        }

      });

    });

    tg.wait();

    return;
}


/* This stage records the number of points per leaf node to create
the CPU and GPU histograms separately */
template<typename pointer_t>
void stageRateStats(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, std::vector<int>& leafs_cpu, std::vector<int>& leafs_gpu, pointer_t qtreeIn, Vector2D min, float rate)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    Vector2D cellCenter;

    int chunk = static_cast<int>(nRows * rate);

    // std::vector<int> leafs_gpu;

    for(int jj = 0 ; jj < chunk ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < nCols ; ii++ ){

            cellCenter.x = initX + ii*Displace;

            register_data(cellCenter, qtreeIn, Wsize*0.5, leafs_gpu);

        }

    }

    // std::vector<int> leafs_cpu;

    for(int jj = chunk ; jj < nRows ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < nCols ; ii++ ){

            cellCenter.x = initX + ii*Displace;

            register_data(cellCenter, qtreeIn, Wsize*0.5, leafs_cpu);

        }

    }

    return;
}


/* sequential stage1 */
template<typename pointer_t>
void stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, pointer_t qtree, Vector2D min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    Vector2D cellCenter;

    Lpoint newmin;

    int cellPoints;

    // Tengo que hacerlo así porque el kernel no captura las variables globales
    // QtreeG5 p_arbol = array_indexes;
    // Lpoint* p_puntos = array_all_points;

    for(int jj = 0 ; jj < Ccol ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < Crow ; ii++ ){

          // printf("CELDA %d,%d <-------------------------------------------------------\n",jj,ii);

            cellCenter.x = initX + ii*Displace;

            // newmin = searchNeighborsMin(&cellCenter, qtree, Wsize/2, &cellPoints);
            // newmin = cpuSearchNeighborsMin(cellCenter, qtree, Wsize/2, cellPoints);
            newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
            // newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

            if(cellPoints >= minNumPoints ){
              // printf("minimo: %d", newmin.id);
              minIDs[jj*Crow+ii] = newmin.id;
            }
        }

    }

    return;
}


/* Stage1 with parallel_reduce */
template<typename pointer_t>
std::vector<int> stage1reduce(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
                        unsigned short minNumPoints, pointer_t qtree, Vector2D min)
{
  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  std::vector<int> init;

  // Tengo que hacerlo así para gpuSearchIndex porque el kernel no captura las variables globales
  // QtreeG5 p_arbol = array_indexes;
  // Lpoint* p_puntos = array_all_points;

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

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                cellCenter.x = initX + ii*Displace;

                // newmin = searchNeighborsMin(cellCenter, qtree, Wsize/2, cellPoints);
                // newmin = searchUpdated(cellCenter, qtree, Wsize*0.5, cellPoints);
                newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
                // Lpoint newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

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


/* Stage1 with parallel_for that remembers previous minimum */
template<typename pointer_t>
void stage1tbbRem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, pointer_t qtreeIn, Vector2D min)
{

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
        // int cellPoints = 0;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            int cellPoints = 0;

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

              cellCenter.x = initX + ii*Displace;

              makeBox(cellCenter, Wsize*0.5, boxMin, boxMax);

              if(insideBox2D(&newmin,boxMin,boxMax)){
              // if(cellPoints > 0 && insideBox2D(&newmin,boxMin,boxMax)){

                Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

                int old_cellPoints = cellPoints;

                Lpoint tmp = gpuSearchNeighborsMin(oCell, qtreeIn, Displace*0.5, Wsize*0.5, cellPoints);

                // We're assuming the points were equidistant throughout the cell, which isn't always true.

                /*En este punto, si queremos ser estrictos, en vez de hacer esta suposición podemos
                lanzar una búsqueda "countMin" en la zona en la que conocemos el mínimo, pudiendo lanzar
                las dos búsquedas diferentes en paralelo con un "parallel_invoke" */
                cellPoints += (int)(old_cellPoints * Overlap);

                if(tmp.z < newmin.z){
                  newmin = tmp;
                }


              } else {

                // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
                newmin = gpuSearchNeighborsMin(cellCenter, qtreeIn, Wsize*0.5, Wsize*0.5, cellPoints);

              }

              if( cellPoints >= minNumPoints ){

                  // v.push_back(newmin.id);
                  minIDs[jj*Crow + ii] = newmin.id;

              }

            }
        }

    });

    return;
}

/* Sequential Stage2 */
uint32_t stage2(uint32_t countMin, std::vector<int>& minIDs){

    uint32_t index = 0;

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


/* Sequential Stage3 */
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

        countNeighbors2D(cellCenter, grid, Bsize/2, cellPoints);

        // Tengo que hacerlo porque el algoritmo no lo hace
        if(cellPoints == 0){

          newmin = searchNeighborsMin(cellCenter, qtreeIn, Bsize/2, cellPoints);

          if(cellPoints>0){
            minGridIDs.push_back(newmin.id);
          }
        }
      }
    }

    return;
}


/* Stage3 with parallel_reduce */
template<typename pointer_t>
std::vector<int> stage3reduce(unsigned short Bsize, unsigned short nCols, unsigned short nRows,
                            pointer_t qtree, Qtree grid, Vector2D min)
{
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

                countNeighbors2D(cellCenter, grid, Bsize*0.5, cellPoints);

                if( cellPoints == 0 ){

#ifdef INDEX
                    newmin = searchNeighborsMin(cellCenter, qtree, Bsize*0.5, cellPoints);
#else                    
                    newmin = gpuSearchNeighborsMin(cellCenter, qtree, Bsize*0.5, Bsize*0.5, cellPoints);
#endif                    
                    

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


/* Saves results in .xyz format */
int save_file(std::string file_name, std::vector<int>& ids, Lpoint** pointer)
{
  // uint32_t size = pointList.size();
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


/* Checks results by comparing with GOLD result */
double check_results(std::string filename, std::vector<int>& ids, Lpoint** pointer, float Displace)
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
    makeBox(aux, Displace*0.5, boxMin, boxMax);
    for(auto pg : vgold){
      if(insideBox2D(pg,boxMin,boxMax)){
        if(fabs(p.z - pg.z) < 0.01){
          count++;
          break;
          }
        }
    }
  }

  double rate = count/((double)(npoints))*100.0;
  printf("%d points correct; %.2f%%\n", count, rate);

  return rate;

  // return std::equal(v.begin(), v.end(), v2.begin(), [](const Vector3D a, const Vector3D b){
  //   return fabs(a.z - b.z) < 0.01;
  // });
}


/* Reads the point cloud in .xyz format */
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


/* Reads the point cloud in .xyz format */
int readXYZfile(std::string filename, Lpoint* & point_cloud, uint64_t & Npoints, Vector2D &min, Vector2D &max) 
{
  FILE* fileXYZ;
  if((fileXYZ = fopen(filename.c_str(),"r")) == NULL){
    printf("Unable to open file!\n");
    return -1;
  }
  if ( filename.find("ArzuaH.xyz") != std::string::npos || filename.find("AlcoyH.xyz") != std::string::npos || 
       filename.find("BrionFH.xyz") != std::string::npos || filename.find("BrionUH.xyz") != std::string::npos ){
    printf("Read header...\n");
    if(fscanf(fileXYZ, "%lu\n%lf\n%lf\n%lf\n%lf\n",&Npoints, &min.x, &max.x, &min.y, &max.y) < 5){
        printf("Imposible to read header values\n");
        return -1;
    }
  }
  // Allocate memory for the LiDAR points
  try {
    point_cloud = static_cast<Lpoint*>(mallocWrap(Npoints * sizeof(Lpoint)));
  } catch (sycl::exception &E) {
    std::cout << "point_cloud malloc: " << E.what() << std::endl;
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

int read_pointsC(std::string file_name, Lpoint* point_cloud)
{
  FILE* fileLAS;

  if((fileLAS = fopen(file_name.c_str(),"r")) == NULL){
    printf("Unable to open file!\n");
    return -1;
  }

  for(int i=0; i<Npoints ; i++){
    //Obtengo los datos id X Y Z
    point_cloud[i].id = i;
    if(fscanf(fileLAS, "%lf %lf %lf",&point_cloud[i].x,&point_cloud[i].y,&point_cloud[i].z) < 3){
      printf("Imposible to obtain values\n");
      return -1;
    }
    while(fgetc(fileLAS)!='\n');
  }

  //Ya no necesito mas el fichero
  if(fclose(fileLAS)){
    printf("Cannot close the file\n");
    return -1;
  }

  return 0;
}


void* mallocWrap(size_t size)
{
  // void *ptr = std::malloc(size);
  void *ptr = malloc_host(size, gpu_queue);
  // void *ptr = sycl::malloc(size, qdevice, qcontext, global_alloc);
  if (ptr)
      return ptr;
  else
      throw std::bad_alloc{};
}


template<typename pointer_t>
void freeWrap(pointer_t& ptr)
{
  // free(ptr);
  free(ptr, gpu_queue);
  ptr = NULL;

  return;
}


/* Recursively copy the tree using pointers */
template<typename cpunode_t>
void copyQtree2(cpunode_t qt, QtreeG4 parent)
{
  if(isLeaf(qt)) {

    // si este nodo es hoja, solo tengo que copiar los puntos ya que lo he guardado anteriormente
    // quadrants lo dejo a NULL porque es hoja; evito ocupar más espacio

    size_t vsize = qt->points.size();

    if(vsize > 0) {

      uint64_t index = num_pts_copiados;

      parent->points = &(array_all_points[num_pts_copiados]);
      
      parent->numPts = vsize;

      num_pts_copiados += vsize;

      Lpoint* min = NULL;
      for(Lpoint* p : qt->points) {
        array_all_points[index] = *p;
      
        if(min == NULL) {
          // min = p;
          min = &(array_all_points[index]);
        } else if(p->z < min->z) {
          min = &(array_all_points[index]);
        }

        index++;
      }

      parent->min = min;
    }

  } else {

    // ahora es cuando se que quadrants debe apuntar a algo
    parent->quadrants = &(array_pointers[copied_nodes]);

    // copio los 4 cuadrantes contiguamente
    for(int i=0; i<4; ++i) {

      cpunode_t src = qt->quadrants[i];

      QtreeG4 dst = &(parent->quadrants[i]);

      new (dst) QtreeG4_t(  parent, 
                            src->center, 
                            src->radius );
      
      // copied_nodes++;
    }

    copied_nodes += 4;

    // continúo recursivamente para cada uno de los cuadrantes
    Lpoint* min = NULL;
    for(int i=0; i<4; ++i) {

      QtreeG4 dst = &(parent->quadrants[i]);

      copyQtree2(qt->quadrants[i], dst);

      parent->numPts += dst->numPts;

      if(dst->min != NULL) {
        if(min == NULL) {
          min = dst->min;
        } else if(dst->min->z < min->z) {
          min = dst->min;
        }
      }

    }

    parent->min = min;

  }

}

template<typename cpunode_t>
QtreeG4 launchCopy(cpunode_t qt)
{
  QtreeG4 root = &(array_pointers[0]);

  new (root) QtreeG4_t(  NULL, 
                        qt->center, 
                        qt->radius );

  copied_nodes++; // hasta aquí ya tengo copiado root

  copyQtree2(qt, root);
  
  return root;
  
}



/* Recursively copy the tree using indexes */
template<typename cpunode_t>
void cudaCopyQtree(cpunode_t qt, int32_t parentID)
{
  QtreeG5 parent = &(array_indexes[parentID]);
  
  if(isLeaf(qt)) {

    // si este nodo es hoja, solo tengo que copiar los puntos ya que lo he guardado anteriormente
    // quadrants lo dejo a NULL porque es hoja

    size_t vsize = qt->points.size();
    if(vsize > 0) {

      parent->points = cuda_num_pts_copiados;

      parent->numPts = vsize;

      int32_t min = -1;
      for(Lpoint* p : qt->points) {

        array_all_points[cuda_num_pts_copiados] = *p;
      
        if(min < 0) {
          min = cuda_num_pts_copiados;
        } else if(p->z < array_all_points[min].z) {
          min = cuda_num_pts_copiados;
        }

        cuda_num_pts_copiados++;
      }
      parent->min = min;
    }

  } else {

    int32_t quadrantsID = cuda_copied_nodes;
    // ahora es cuando se que quadrants debe apuntar a algo
    parent->quadrants = cuda_copied_nodes;

    // copio los 4 cuadrantes contiguamente
    for(int i=0; i<4; ++i) {

      cpunode_t src = qt->quadrants[i];

      QtreeG5 dst = &(array_indexes[quadrantsID + i]);

      new (dst) QtreeG5_t(  parentID, 
                            src->center, 
                            src->radius );
      
      cuda_copied_nodes++;
    }

    // continúo recursivamente para cada uno de los cuadrantes
    int32_t min = -1;
    for(int i=0; i<4; ++i) {

      QtreeG5 dst = &(array_indexes[quadrantsID + i]);

      cudaCopyQtree(qt->quadrants[i], quadrantsID + i);

      parent->numPts += dst->numPts;

      if(dst->min > -1) {
        if(min < 0) {
          min = dst->min;
        } else if(array_all_points[dst->min].z < array_all_points[min].z) {
          min = dst->min;
        }
      }

    }
    parent->min = min;

  }

}

template<typename cpunode_t>
void cudaQtree(cpunode_t qt)
{
  QtreeG5 root = &(array_indexes[0]);

  new (root) QtreeG5_t( 0, 
                        qt->center, 
                        qt->radius );

  cuda_copied_nodes++; // hasta aquí ya tengo copiado root

  cudaCopyQtree(qt, 0);
  
  return;
}



/* Con este método de copia lo que hago es alinear los nodos que voy a ir
consultando conforme avanza el método de búsqueda de forma natural;
El problema de este método es que, al no estar los nodos contiguos en memoria,
ya no me vale con mantener en el nodo guardada solo la posición del primero,
sino que necesito guardar la posición de los cuatro, por tanto, el manejo 
de la reserva de memoria con USM debe ser diferente */

// template<typename cpunode_t>
// QtreeG4 copyQtreeP(cpunode_t qt, QtreeG4 parent)
// {
//   QtreeG4 pos;

//   if(isLeaf(qt)) {

//     /*Si este nodo de CPU es hoja, tengo que guardar el nodo y copiar los puntos*/
//     pos = &(array_pointers[copied_nodes]);

//     new (pos) QtreeG4_t(  parent, 
//                           qt->center, 
//                           qt->radius, 
//                           p_alloc );

//     copied_nodes++;

//     size_t vsize = qt->points.size();

//     if(vsize > 0) {
//       pos->points = &(array_all_points[num_pts_copiados]);
//       // size_t vsize = qt->points.size();
//       pos->numPts = vsize;

//       // std::memcpy(&(array_all_points[num_pts_copiados]), qt->points.data(), pos->numPts*sizeof(Lpoint*));
//       // num_pts_copiados += pos->numPts;

//       Lpoint* min = NULL;
//       for(Lpoint* p : qt->points) {
//         array_all_points[num_pts_copiados] = *p;
//         // num_pts_copiados++;
      
//         // Lpoint** min = &(pos->min);
//         if(min == NULL) {
//           // min = p;
//           min = &(array_all_points[num_pts_copiados]);
//         } else if(p->z < min->z) {
//           min = &(array_all_points[num_pts_copiados]);
//         }

//         num_pts_copiados++;
//       }
//       pos->min = min;
//     }

//   } else {

//     pos = &(array_pointers[copied_nodes]);

//     new (pos) QtreeG4_t(  parent, 
//                           qt->center, 
//                           qt->radius, 
//                           p_alloc );

//     copied_nodes++;

//     /* aqui falta actualizar las posiciones de los cuadrantes
//     en este caso lo cuadrantes NO van a estar contiguos en memoria*/
//     Lpoint* min = NULL;
//     for(int i=0; i<4; ++i) {

//       QtreeG4 newqt = copyQtreeP(qt->quadrants[i], pos);

//       pos->quadrants[i] = newqt;

//       pos->numPts += newqt->numPts;

//       if(newqt->min != NULL) {
//         if(min == NULL) {
//           min = newqt->min;
//         } else if(newqt->min->z < min->z) {
//           min = newqt->min;
//         }
//       }

//     }
//     pos->min = min;

//   }

//   return pos;

// }


/*Este se utiliza cuando NO tengo los cuadrantes consecuivos en memoria*/
// template<typename pointer_t>
// void findPosition(pointer_t parent, pointer_t& current, int& idx) 
// {

//   if(parent->quadrants[0] == current){
//     // std::advance(iterador, 1);
//     idx=1;
//     current = parent->quadrants[1];
//   }
//   else if(parent->quadrants[1] == current){
//     // std::advance(iterador, 2);
//     idx=2;
//     current = parent->quadrants[2];
//   }
//   else if(parent->quadrants[2] == current){
//     // std::advance(iterador, 3);
//     idx=3;
//     current = parent->quadrants[3];
//   } 
//   else {
//     //lo pongo apuntando al último
//     // std::advance(iterador, 3);
//     idx = 4; // termino, pero dejo listo current
//     current = parent->quadrants[3];
//   }

// }

/*Este se utiliza cuando NO tengo los cuadrantes consecuivos en memoria*/
// Lpoint gpuSearchNeighborsMin(Vector2D& center, QtreeG4 qtree, float radius, int& numNeighs)
// {
//     int idx = 0;
//     int numInside = 0;
//     Vector2D boxMin, boxMax;

//     // *numNeighs = 0;
//     makeBox(center, radius, boxMin, boxMax);

//     Lpoint min = {0,0.0,0.0,std::numeric_limits<double>::max()};

//     QtreeG4 parent = qtree;
//     QtreeG4 current = qtree->quadrants[0];

//     while(current != qtree) {


//       if(idx > 3) {

//         current = current->parent;

//         if(current != qtree){

//           parent = current->parent;

//           current = findPosition(parent, current, idx);

//         }
        
//       } else {

//         if( current->quadrants[0] == NULL ) { // isLeaf??

//           if( current->numPts > 0 ) { //isEmpty??

//             /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
//             if(boxInside2D(boxMin, boxMax, current)) {

//               if(current->min->z < min.z)
//                 min = *(current->min);
//               numInside += current->numPts;

//             /*Si está parcialmente solapada, tengo dos opciones*/
//             } else if(boxOverlap2D(boxMin, boxMax, current)) {

//               int N = current->numPts;
//               Lpoint* pointer = current->points;

//               /*Si el mínimo está dentro solo comparo una vez y cuento */
//               if( insideBox2D(current->min, boxMin, boxMax) ) {

//                 if(current->min->z < min.z)
//                   min = *(current->min);

//                 for(int i=0; i<N; ++i) {
//                   if(insideBox2D(&pointer[i], boxMin, boxMax))
//                   {
//                     numInside++;
//                   }
//                 }
//               }
//               /*Si el mínimo está fuera, tengo que comprobarlo todo*/
//               else {
//                 for(int i=0; i<N; ++i) {
//                   if(insideBox2D(&pointer[i], boxMin, boxMax))
//                   {
//                     if (pointer[i].z < min.z) {
//                         min = pointer[i];
//                     }
//                     numInside++;
//                   }
//                 }
//               }

//             }

//           }

//           idx++;
//           if(idx < 4) {
//             current = parent->quadrants[idx];

//           }

//         } else {
          
//           if(!boxOverlap2D(boxMin, boxMax, current) || current->min == NULL) { //No solapada o vacia

//             idx++;
//             if(idx < 4) {
//               current = parent->quadrants[idx];

//             }
//           }
//           /*si la caja está completamente solapada, capturo el mínmo y el numero de puntos, y continuo*/
//           else if( boxInside2D(boxMin, boxMax, current) ) {

//             if(current->min->z < min.z)
//               min = *(current->min);
//             numInside += current->numPts;

//             idx++;
//             if(idx < 4) {
//               current = parent->quadrants[idx];

//             }
//           }
//           /*pero si solo está parcialmente solapada, profundizo*/
//           else {

//             idx = 0;
//             parent = current;
//             current = current->quadrants[0];

//           }
//         }
//       }

//     }

//     numNeighs = numInside;

//     return min;
// }


/* Performs the search when the quadrants are contiguos in memory */
Lpoint gpuSearchNeighborsMin(Vector2D& center, QtreeG4 qtree, float radiusX, float radiusY, int& numNeighs)
{
    int idx = 0;
    int numInside = 0;
    Vector2D boxMin, boxMax;

    // *numNeighs = 0;
    makeBox(center, radiusX, radiusY, boxMin, boxMax);

    Lpoint min = {0,0.0,0.0,std::numeric_limits<double>::max()};

    QtreeG4 parent = qtree;
    QtreeG4 current = qtree->quadrants;

    while(current != qtree) {


      if(idx > 3) {

        current = current->parent;

        if(current != qtree){

          parent = current->parent;

          findPosition(parent, current, idx);

        }
        
      } else {

        if( current->quadrants == NULL ) { // isLeaf??

          if( current->numPts > 0 ) { //isEmpty??

            /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
            if(boxInside2D(boxMin, boxMax, current)) {

              if(current->min->z < min.z)
                min = *(current->min);
              numInside += current->numPts;

            /*Si está parcialmente solapada, tengo dos opciones*/
            } else if(boxOverlap2D(boxMin, boxMax, current)) {

              int N = current->numPts;
              Lpoint* pointer = current->points;

              /*Si el mínimo está dentro solo comparo una vez y cuento */
              if( insideBox2D(current->min, boxMin, boxMax) ) {

                if(current->min->z < min.z)
                  min = *(current->min);

                for(int i=0; i<N; ++i) {
                  if(insideBox2D(&pointer[i], boxMin, boxMax))
                  {
                    numInside++;
                  }
                }
              }
              /*Si el mínimo está fuera, tengo que comprobarlo todo*/
              else {
                for(int i=0; i<N; ++i) {
                  if(insideBox2D(&pointer[i], boxMin, boxMax))
                  {
                    if (pointer[i].z < min.z) {
                        min = pointer[i];
                    }
                    numInside++;
                  }
                }
              }

            }
            // else {
            //   std::cout << "Ni lo uno ni lo otro\n" << std::flush;
            // }

          }

          idx++;
          if(idx < 4) {
            current = &(parent->quadrants[idx]);

          }

        } else {
          
          if(!boxOverlap2D(boxMin, boxMax, current) || current->numPts == 0) { //No solapada o vacia

            idx++;
            if(idx < 4) {
              current = &(parent->quadrants[idx]);

            }
          }
          /*si la caja está completamente solapada, capturo el mínmo
           y el numero de puntos, y continuo*/
          else if( boxInside2D(boxMin, boxMax, current) ) {

            if(current->min->z < min.z)
              min = *(current->min);
            numInside += current->numPts;

            idx++;
            if(idx < 4) {
              current = &(parent->quadrants[idx]);

            }
          }
          /*pero si solo está parcialmente solapada, profundizo*/
          else {

            idx = 0;
            parent = current;
            current = current->quadrants;

          }
        }
      }

    }

    numNeighs = numInside;

    return min;
}

template<typename pointer_t>
Lpoint findUpdated(pointer_t qtree, Vector2D& boxMin, Vector2D& boxMax, int& numInside)
{
    // Lpoint tmp, min = nomin;
    Lpoint tmp, min = {0,0.0,0.0,std::numeric_limits<double>::max()};

    if(qtree->quadrants == NULL)
    {

      // if(qtree->numPts == 0)
      //   std::cout << "0 puntos\n" << std::flush;

      int N = qtree->numPts;
      Lpoint* pointer = qtree->points;

      /*Si el mínimo está dentro solo comparo una vez y cuento */
      if( insideBox2D(qtree->min, boxMin, boxMax) ) {

        if(qtree->min->z < min.z)
          min = *(qtree->min);

        for(int i=0; i<N; ++i) {
          if(insideBox2D(&pointer[i], boxMin, boxMax))
          {
            numInside++;
          }
        }
      }
      /*Si el mínimo está fuera, tengo que comprobarlo todo*/
      else {
        for(int i=0; i<N; ++i) {
          if(insideBox2D(&pointer[i], boxMin, boxMax))
          {
            if (pointer[i].z < min.z) {
                min = pointer[i];
            }
            numInside++;
          }
        }
      }

    } else {


      for(int i = 0; i < 4; i++) {

        pointer_t current = &(qtree->quadrants[i]);
        // Check
        if(!boxOverlap2D(boxMin, boxMax, current) || current->numPts == 0)
            continue;

        // completely overlaped
        else if( boxInside2D(boxMin, boxMax, current) ) {

          if(current->min->z < min.z)
            min = *(current->min);
          numInside += current->numPts;

        }
        // At this point it's not empty and it's overlapped
        else {
          tmp = findUpdated(current, boxMin, boxMax, numInside);
          if (tmp.z < min.z) {
              min = tmp;
          }
        }

      } //for

    }

    return min;
}

template<typename pointer_t>
Lpoint searchUpdated(Vector2D& center, pointer_t qtree, float radius, int& numNeighs)
{
    Vector2D boxMin, boxMax;

    numNeighs = 0;
    makeBox(center, radius, boxMin, boxMax);

    return findUpdated(qtree, boxMin, boxMax, numNeighs);
}


/* It looks for the posisition of the current node in the 
quadrant array of the parent node */
template<typename accNode_t>
void findPosition(accNode_t& arrayNodes, QtreeG5 parent, QtreeG5& current, int& idx) 
{

  QtreeG5 quadrants = &(arrayNodes[parent->quadrants]);

  if(quadrants == current){
    // std::advance(iterador, 1);
    idx=1;
    current = &(quadrants[1]);
  }
  else if(&(quadrants[1]) == current){
    // std::advance(iterador, 2);
    idx=2;
    current = &(quadrants[2]);
  }
  else if(&(quadrants[2]) == current){
    // std::advance(iterador, 3);
    idx=3;
    current = &(quadrants[3]);
  } 
  else {
    //lo pongo apuntando al último
    // std::advance(iterador, 3);
    idx = 4; // termino, pero dejo listo current
    current = &(quadrants[3]);
  }

}

/* Performs the search when the quadrants are contiguos in memory
and indexes are used */
template<typename accNode_t, typename accPoint_t>
Lpoint gpuSearchIndex( Vector2D& center, accNode_t aNodes, accPoint_t aPoints, float radiusX, float radiusY, int& numNeighs)
{
    int idx = 0;
    int numInside = 0;
    Vector2D boxMin, boxMax;

    makeBox(center, radiusX, radiusY, boxMin, boxMax);

    Lpoint min = {0,0.0,0.0,std::numeric_limits<double>::max()};

    QtreeG5 root = &(aNodes[0]);
    
    QtreeG5 parent = root;
    QtreeG5 current = &(aNodes[root->quadrants]);

    while(current != root) {

      if(idx > 3) {

        current = &(aNodes[current->parent]);

        if(current != root){

          parent = &(aNodes[current->parent]);

          // current = findPosition(aNodes, parent, current, idx);
          findPosition(aNodes, parent, current, idx);

        }
        
      } else {

        if( current->quadrants < 0 ) { // isLeaf??

          if( current->numPts > 0 ) { //isEmpty??

            /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
            if(boxInside2D(boxMin, boxMax, current)) {

              Lpoint* currentMin = &(aPoints[current->min]);

              if(currentMin->z < min.z)
                min = *currentMin;
              numInside += current->numPts;

            /*Si está parcialmente solapada, tengo dos opciones*/
            } else if(boxOverlap2D(boxMin, boxMax, current)) {

              int N = current->numPts;

              Lpoint* pointer = &(aPoints[current->points]);

              /*Si el mínimo está dentro solo comparo una vez y cuento */
              Lpoint* currentMin = &(aPoints[current->min]);

              if( insideBox2D(currentMin, boxMin, boxMax) ) {

                if(currentMin->z < min.z)
                  min = *currentMin;

                for(int i=0; i<N; ++i) {
                  if(insideBox2D(&pointer[i], boxMin, boxMax))
                  {
                    numInside++;
                  }
                }
              }
              /*Si el mínimo está fuera, tengo que comprobarlo todo*/
              else {
                for(int i=0; i<N; ++i) {
                  if(insideBox2D(&pointer[i], boxMin, boxMax))
                  {
                    if (pointer[i].z < min.z) {
                        min = pointer[i];
                    }
                    numInside++;
                  }
                }
              }

            }

          }

          idx++;
          if(idx < 4) {

            current = &(aNodes[parent->quadrants + idx]);

          }

        } else {
          
          if(!boxOverlap2D(boxMin, boxMax, current) || current->numPts == 0) { //No solapada o vacia

            idx++;
            if(idx < 4) {
              current = &(aNodes[parent->quadrants + idx]);

            }
          }

          else if( boxInside2D(boxMin, boxMax, current) ) {

            Lpoint* currentMin = &(aPoints[current->min]);

            if(currentMin->z < min.z)
              min = *currentMin;
            numInside += current->numPts;

            idx++;
            if(idx < 4) {
              current = &(aNodes[parent->quadrants + idx]);

            }
          }
          /*pero si solo está parcialmente solapada, profundizo*/
          else {

            idx = 0;
            parent = current;
            current = &(aNodes[current->quadrants]);

          }
        }
      }

    }

    numNeighs = numInside;

    return min;
}



/* It only processes the point cloud on GPUs */
template<typename pointer_t>
void stage1gpu(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, pointer_t qtree, Vector2D min)
{
    is_gpu_used = 1;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // Tengo que hacerlo así porque el kernel no captura las variables globales
    // QtreeG5 p_arbol = array_indexes;
    // Lpoint* p_puntos = array_all_points;

    gpu_queue.parallel_for(range<2>(Ccol,Crow), [=](id<2> index) { 
    // gpu_queue.parallel_for(range<1>(Ccol), [=](id<1> j) { 

      int jj = static_cast<int>(index[0]);
      int ii = static_cast<int>(index[1]);

      Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

      int cellPoints;

      Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
      // Lpoint newmin = gpuSearchNeighborsMinSimple(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
      // Lpoint newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

      if( cellPoints >= minNumPoints ){
        minGPU[jj*Crow+ii] = newmin.id;
      }

    }).wait();

    return;
}


template<typename pointer_t>
void stage1gpuRem(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, int* minIDs, pointer_t qtree, Vector2D min)
{

    // std::cout << " Using device: " << gpu_queue.get_device().get_info<info::device::name>() << std::endl;
    is_gpu_used = 1;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    gpu_queue.parallel_for(range<2>(nRows, nCols), [=](sycl::id<2> index) { 

      int jj = static_cast<int>(index[0]);
      int ii = static_cast<int>(index[1]);
      int cellPoints = 0;
      Vector2D boxMax, boxMin;
      Lpoint newmin = {0,0.0,0.0,0.0};

      Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

      makeBox(cellCenter, Wsize*0.5, boxMin, boxMax);

      if(insideBox2D(&newmin,boxMin,boxMax)){ // is the latest point inside the new box?

        Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

        int old_cellPoints = cellPoints;

        Lpoint tmp = gpuSearchNeighborsMin(oCell, qtree, Displace*0.5, Wsize*0.5, cellPoints);

        // We're assuming the points were equidistant throughout the cell, which isn't always true.
        cellPoints += (int)(old_cellPoints * Overlap);

        if(tmp.z < newmin.z){
          newmin = tmp;
        }

      } else {

        newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);

      }

      if( cellPoints >= minNumPoints ){

          // v.push_back(newmin.id);
          minIDs[jj*nCols + ii] = newmin.id;

      }


    }).wait();

    return;
}

void stage1index(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, Vector2D min)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, global_alloc>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    // std::cout << " Using device: " << gpu_queue.get_device().get_info<info::device::name>() << std::endl;

    double Displace = round2d(Wsize*(1.-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    auto nR = range<1>(cuda_copied_nodes);
    auto pR = range<1>(cuda_num_pts_copiados);
    auto mR = range<1>(Ccol*Crow);

    // printf("Vector nodos: %d, y de puntos: %d", cuda_copied_nodes, cuda_num_pts_copiados);

    buffer<QtreeG5_t, 1> bufNodes(array_indexes, nR);
    buffer<Lpoint, 1> bufPoints(array_all_points, pR);
    buffer<int, 1> bufMins(minGPU, mR);

    gpu_queue.submit([&](handler &h) {

      /*Hay que poner modo read_write porque de no ser así da un error: 
      " cannot initialize a variable of type 'const QtreeG5' (aka 'QtreeG5_t *const') 
      with an rvalue of type 'const QtreeG5_t *' " cuando se llama a gpuSearchIndex */

      auto accNodes = bufNodes.get_access<access::mode::read_write>(h);
      auto accPoints = bufPoints.get_access<access::mode::read_write>(h);
      auto accC = bufMins.get_access<access::mode::write>(h);

      // sycl::stream out(1024, 256, h);

      // # nd-range kernel
      h.parallel_for(range<2>(Ccol,Crow), [=](id<2> index) { 
      // h.parallel_for(range<1>(Ccol), [=](id<1> j) { 

        int jj = static_cast<int>(index[0]);
        int ii = static_cast<int>(index[1]);

        Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

        int cellPoints;

        Lpoint newmin = gpuSearchIndex(cellCenter, accNodes, accPoints, Wsize*0.5, Wsize*0.5, cellPoints);

        // out << " (" << Ccol << ", " << Crow << ") " <<  "minimo: (" << ii << ", " << jj << ") " << endl;

        if( cellPoints >= minNumPoints ){

            // out << "minimo: (" << ii << ", " << jj << ") " << endl;
            // out << "minimo: " << newmin.id << endl;
            // minGPU[jj*Crow+ii] = newmin.id;
            accC[jj*Crow+ii] = newmin.id;

            // v.push_back(newmin.id);
            // printf("Tengo un minimo!\n");
            // std::cout << "Tengo un minimo\n";

        }

      });

    // }).wait();
    });
    // bloqueante; ya deberíane estar los datos en minGPU
    auto hMins = bufMins.get_access<access::mode::read>();

    // std::memcpy(minGPU, bufMins, Ccol*Crow*sizeof(int));

    return;
}

void stage1heterAcc(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, int* minGPU, Vector2D min, float rate)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, global_alloc>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    // std::cout << " Using device: " << gpu_queue.get_device().get_info<info::device::name>() << std::endl;

    is_gpu_used = 1;

    double Displace = round2d(Wsize*(1.-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    int chunk = static_cast<int>(nRows * rate);

    auto nR = range<1>(cuda_copied_nodes);
    auto pR = range<1>(cuda_num_pts_copiados);
    // auto mR = range<1>(nRows*nCols);
    auto mR = range<1>(chunk*nCols);

    // printf("Vector nodos: %d, y de puntos: %d", cuda_copied_nodes, cuda_num_pts_copiados);

    buffer<QtreeG5_t, 1> bufNodes(array_indexes, nR);
    buffer<Lpoint, 1> bufPoints(array_all_points, pR);
    buffer<int, 1> bufMins(minGPU, mR);

    tbb::task_group tg;

    tg.run([&]{

      gpu_queue.submit([&](handler &h) {

        /*Hay que poner modo read_write porque de no ser así da un error: 
        " cannot initialize a variable of type 'const QtreeG5' (aka 'QtreeG5_t *const') 
        with an rvalue of type 'const QtreeG5_t *' " cuando se llama a gpuSearchIndex */

        auto accNodes = bufNodes.get_access<access::mode::read_write>(h);
        auto accPoints = bufPoints.get_access<access::mode::read_write>(h);
        auto accC = bufMins.get_access<access::mode::write>(h);

        // sycl::stream out(1024, 256, h);

        // # nd-range kernel
        h.parallel_for(range<2>(chunk,nCols), [=](id<2> index) { 
        // h.parallel_for(range<1>(nRows), [=](id<1> j) { 

          int jj = static_cast<int>(index[0]);
          int ii = static_cast<int>(index[1]);

          Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

          int cellPoints;

          Lpoint newmin = gpuSearchIndex(cellCenter, accNodes, accPoints, Wsize*0.5, Wsize*0.5, cellPoints);

          // out << " (" << nRows << ", " << nCols << ") " <<  "minimo: (" << ii << ", " << jj << ") " << endl;

          if( cellPoints >= minNumPoints ){

              // out << "minimo: (" << ii << ", " << jj << ") " << endl;
              // out << "minimo: " << newmin.id << endl;
              // minGPU[jj*nCols+ii] = newmin.id;
              accC[jj*nCols+ii] = newmin.id;

              // v.push_back(newmin.id);
              // printf("Tengo un minimo!\n");
              // std::cout << "Tengo un minimo\n";

          }

        });

      // }).wait();
      });

    });

    tg.run([&]{

      // Tengo que hacerlo así porque el kernel no captura las variables globales
      QtreeG5 p_arbol = array_indexes;
      Lpoint* p_puntos = array_all_points;

      tbb::parallel_for( tbb::blocked_range2d<int,int>{chunk,nRows,0,nCols},
                        [&](tbb::blocked_range2d<int,int> r ) {

          // Lpoint newmin;
          Vector2D cellCenter;
          int cellPoints;
          int je = r.rows().end();
          int ie = r.cols().end();

          for(int jj = r.rows().begin(); jj < je; ++jj) {

              cellCenter.y = initY + jj*Displace;

              for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                  cellCenter.x = initX + ii*Displace;

                  // Lpoint newmin = cpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);
                  // Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
                  // newmin = searchNeighborsMin(cellCenter, cpuqt, Wsize/2, cellPoints);
                  // Lpoint newmin = searchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);
                  Lpoint newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

                  if( cellPoints >= minNumPoints ){

                      minGPU[jj*nCols+ii] = newmin.id;

                  }
              }
          }

      });

      // bloqueante; ya deberíane estar los datos en minGPU
      auto hMins = bufMins.get_access<access::mode::read>();

    });

    tg.wait();

    return;
}


template<typename pointer_t>
void stage1heter(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, int* minGPU, pointer_t qtree, Vector2D min, float rate)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, global_alloc>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    // std::cout << " Using device: " << qdevice.get_info<info::device::name>() << std::endl;
    // std::cout << " Using device: " << hqdevice.get_info<info::device::name>() << std::endl;
    is_gpu_used = 1;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // std::vector<int> values(nRows);
    // std::iota(std::begin(values), std::end(values), 0);

    /*Porcentaje de trabajo para GPU*/
    // double rate = 0.2;
    int chunk = static_cast<int>(nRows * rate);

    // gpu_queue.prefetch(array_pointers, cpu_tree_nodes * sizeof(QtreeG4_t));
    // gpu_queue.prefetch(array_all_points, Npoints * sizeof(Lpoint));

    // Tengo que hacerlo así porque el kernel no captura las variables globales
    // QtreeG5 p_arbol = array_indexes;
    // Lpoint* p_puntos = array_all_points;

    sycl::event e = gpu_queue.parallel_for(range<2>(chunk, nCols), [=](sycl::id<2> index) { 
    // gpu_queue.parallel_for(range<1>(nRows), [=](id<1> j) { 

      int jj = static_cast<int>(index[0]);
      int ii = static_cast<int>(index[1]);

      Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

      int cellPoints;

      Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
      // Lpoint newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);


      if( cellPoints >= minNumPoints ){
        minGPU[jj*nCols+ii] = newmin.id;
      }

    });

    tbb::parallel_for( tbb::blocked_range2d<int,int>{chunk,nRows,0,nCols},
                      [&](tbb::blocked_range2d<int,int> r ) {

        // Lpoint newmin;
        Vector2D cellCenter;
        int cellPoints;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                cellCenter.x = initX + ii*Displace;

                // Lpoint newmin = cpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);
                Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
                // newmin = searchNeighborsMin(cellCenter, cpuqt, Wsize/2, cellPoints);
                // Lpoint newmin = searchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);
                // Lpoint newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

                if( cellPoints >= minNumPoints ){

                    minGPU[jj*nCols+ii] = newmin.id;

                }
            }
        }

    });

    e.wait();

    return;
}

template<typename pointer_t>
void stage1hetertg(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, int* minGPU, pointer_t qtree, Vector2D min, float rate)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, global_alloc>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    // std::cout << " Using device: " << qdevice.get_info<info::device::name>() << std::endl;
    // std::cout << " Using device: " << hqdevice.get_info<info::device::name>() << std::endl;
    is_gpu_used = 1;

    double Displace = round2d(Wsize*(1.-Overlap));

    double initX = min.x - Wsize*0.5 + Displace;
    double initY = min.y - Wsize*0.5 + Displace;

    // std::vector<int> values(nRows);
    // std::iota(std::begin(values), std::end(values), 0);

    // tbb::task_group tg;

    /*Porcentaje de trabajo para GPU*/
    // double rate = 0.2;
    int chunk = static_cast<int>(nRows * rate);

    // gpu_queue.prefetch(array_pointers, cpu_tree_nodes * sizeof(QtreeG4_t));
    // gpu_queue.prefetch(array_all_points, Npoints * sizeof(Lpoint));

    // Tengo que hacerlo así porque el kernel no captura las variables globales
    // QtreeG5 p_arbol = array_indexes;
    // Lpoint* p_puntos = array_all_points;

    global_tg.run( [&]{

      gpu_queue.parallel_for(range<2>(chunk, nCols), [=](sycl::id<2> index) { 
      // gpu_queue.parallel_for(range<1>(nRows), [=](id<1> j) { 

        int jj = static_cast<int>(index[0]);
        int ii = static_cast<int>(index[1]);

        Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

        int cellPoints;

        Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
        // Lpoint newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

        if( cellPoints >= minNumPoints ){
          minGPU[jj*nCols+ii] = newmin.id;
        }

      }).wait();
      // }).wait_and_throw();

    });

    global_tg.run( [&]{

      tbb::parallel_for( tbb::blocked_range2d<int,int>{chunk,nRows,0,nCols},
                        [&](tbb::blocked_range2d<int,int> r ) {

          // Lpoint newmin;
          Vector2D cellCenter;
          int cellPoints;
          int je = r.rows().end();
          int ie = r.cols().end();

          for(int jj = r.rows().begin(); jj < je; ++jj) {

              cellCenter.y = initY + jj*Displace;

              for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                  cellCenter.x = initX + ii*Displace;

                  // Lpoint newmin = cpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);
                  Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
                  // newmin = searchNeighborsMin(cellCenter, cpuqt, Wsize/2, cellPoints);
                  // Lpoint newmin = searchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);
                  // Lpoint newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

                  if( cellPoints >= minNumPoints ){

                      minGPU[jj*nCols+ii] = newmin.id;

                  }
              }
          }

      });

    });

    global_tg.wait();

    return;
}

void stage1heterTgMix(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Vector2D min, float rate)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // tbb::blocked_range<int>::size_type gs = 4;

    // tbb::affinity_partitioner aff_p;

    int chunk = static_cast<int>(Ccol * rate);

    tbb::task_group tg;

    // Tengo que hacerlo así porque el kernel no captura las variables globales
    QtreeG5 p_arbol = array_indexes;
    Lpoint* p_puntos = array_all_points;

    tg.run([&]{

      gpu_queue.parallel_for(range<2>(chunk, Crow), [=](sycl::id<2> index) { 

        int jj = static_cast<int>(index[0]);
        int ii = static_cast<int>(index[1]);

        Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

        int cellPoints;

        // Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);
        Lpoint newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

        if( cellPoints >= minNumPoints ){
          minIDs[jj*Crow+ii] = newmin.id;
        }

      }).wait();

    });

    tg.run([&]{

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

                makeBox(cellCenter, Wsize*0.5, boxMin, boxMax);

                if(insideBox2D(&newmin,boxMin,boxMax)){

                  Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

                  int old_cellPoints = cellPoints;

                  // Lpoint tmp = gpuSearchNeighborsMin(oCell, qtreeIn, Displace*0.5, Wsize*0.5, cellPoints);
                  Lpoint tmp = gpuSearchIndex(oCell, p_arbol, p_puntos, Displace*0.5, Wsize*0.5, cellPoints);

                  // We're assuming the points were equidistant throughout the cell, which isn't always true.
                  cellPoints += (int)(old_cellPoints * Overlap);

                  if(tmp.z < newmin.z){
                    newmin = tmp;
                  }

                } else {

                  // newmin = gpuSearchNeighborsMin(cellCenter, qtreeIn, Wsize*0.5, cellPoints);
                  newmin = gpuSearchIndex(cellCenter, p_arbol, p_puntos, Wsize*0.5, Wsize*0.5, cellPoints);

                }

                if( cellPoints >= minNumPoints ){

                    // v.push_back(newmin.id);
                    minIDs[jj*Crow + ii] = newmin.id;

                }

              }
          }

      });

    });

    tg.wait();

    return;
}


template<typename pointer_t>
void stage1heter2q(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, int* minGPU, pointer_t qtree, Vector2D min, float rate)
{

    // std::cout << " Using device: " << qdevice.get_info<info::device::name>() << std::endl;
    // std::cout << " Using device: " << hqdevice.get_info<info::device::name>() << std::endl;

    is_gpu_used = 1;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    /*Porcentaje de trabajo para GPU*/
    int chunk = static_cast<int>(nRows * rate);

    // CPU queue
    sycl::queue hq{cpu_selector_v};

    // gpu_queue.prefetch(array_pointers, cpu_tree_nodes * sizeof(QtreeG4_t));
    // gpu_queue.prefetch(array_all_points, Npoints * sizeof(Lpoint));

    sycl::event e1 = gpu_queue.parallel_for(range<2>(chunk, nCols), [=](sycl::id<2> index) { 
    // sycl::event e1 = gpu_queue.parallel_for(range<2>(chunk, nCols), [=](sycl::id<2> index) { 

      int jj = static_cast<int>(index[0]);
      int ii = static_cast<int>(index[1]);

      Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

      int cellPoints;

      Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);


      if( cellPoints >= minNumPoints ){
        minGPU[jj*nCols+ii] = newmin.id;
      }

    });

    sycl::event e2 = hq.parallel_for(range<2>(nRows-chunk, nCols), [=](sycl::id<2> index) { 
    // sycl::event e2 = hq.parallel_for<class hminSearch>(range<2>(nRows-chunk, nCols), [=](sycl::id<2> index) { 

      int jj = static_cast<int>(index[0]) + chunk;
      int ii = static_cast<int>(index[1]);

      Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

      int cellPoints;

      Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);


      if( cellPoints >= minNumPoints ){
        minGPU[jj*nCols+ii] = newmin.id;
      }

    });

    e1.wait();
    e2.wait();

    return;
}


uint32_t stage2GPU(uint32_t countMin, int* minIDs){

  uint32_t index = 0;

  int ii,jj,id;

  /* esta es la forma de descartar las posiciones no utilizadas, inicializadas a -1 */
  for( jj=0 ;  minIDs[jj] == -1 ; ++jj );

  for( ii=jj ; ii<countMin ; ii=jj ){

    id = minIDs[ii];

    for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

    if(jj-ii > 1){
        // if(jj-ii > 1){
        minIDs[index]=id;
        index++;
        // }
    }
  }

    return index;
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
  int ii, jj, id, reps;
  size_t endCount = pCount.size();

  for( int ii=jj ; ii<endCount ; ii=jj ) {

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
  int ii, jj, id, reps;
  size_t endCount = v.size();

  printf("Se han analizado %zu nodos hoja\n", endCount);

  for( int ii=jj ; ii<endCount ; ii=jj ) {

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


/* Saves the execution time and other markers of the algorithm  */
int save_time(std::string file_name, std::string map_name, int numthreads, 
  float minRadius, int maxNumber, int level, double tree_time, double copytree_time, 
  double owm_time, double GPUratio, uint32_t chunkGPU, double correctness)
{
  // uint32_t size = pointList.size();
  std::fstream out(file_name, out.out | out.app);
  // std::ofstream out;
  // out.precision(std::numeric_limits<double>::digits10);
  // out.open(file_name);
  // Point aux;
  if(!out.is_open())
    return -1;

  out << map_name << " " << numthreads << " ";
  // out.precision(3);
  out << std::defaultfloat << minRadius << " " << maxNumber << " " << level << " ";
  // out.precision(6);
  out << std::fixed << tree_time << " " << copytree_time << " " << owm_time << " ";
  // out.precision(1);
  out << is_gpu_used << " " << chunkGPU << " ";
  out << std::defaultfloat << GPUratio << " " << correctness << std::endl;

  out.close();
  if(out.is_open())
    return -1;
  return 0;
}

#endif // ENVI_GPU_HPP