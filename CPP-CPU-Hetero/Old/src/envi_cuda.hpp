#ifndef ENVI_GPU_HPP

#define ENVI_GPU_HPP

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
// #include <atomic>
#include <fstream>
#include <sstream>

#include <CL/sycl.hpp>
#include <CL/sycl/backend/cuda.hpp>
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for.h"
// #include "tbb/concurrent_vector.h"
#include "tbb/blocked_range2d.h"

using namespace cl::sycl;

// NEW TYPES


typedef struct
{
    double x;
    double y;
    // double z;

} Vector2D;

typedef struct
{
    double x;
    double y;
    double z;

} Vector3D;

typedef struct
{
    unsigned int id;
    double x;
    double y;
    double z;

} Lpoint;


struct Qtree_t;

typedef struct Qtree_t* Qtree;

struct Qtree_t {
  // Qtree quadrants[4];
  std::vector<Qtree> quadrants;
  Qtree parent;
  // Vector3D center;
  Vector2D center;
  float radius;
  std::vector<Lpoint*> points;
};

const usm::alloc global_alloc = usm::alloc::host;
// const usm::alloc global_alloc = usm::alloc::shared;
// const usm::alloc global_alloc = usm::alloc::device;


// GLOBAL DEFINITIONS

class CUDASelector : public sycl::device_selector {
public:
  int operator()(const sycl::device &Device) const override {
    using namespace sycl::info;

    const std::string DriverVersion = Device.get_info<cl::sycl::info::device::driver_version>();

    if (Device.is_gpu() && (DriverVersion.find("CUDA") != std::string::npos)) {
      std::cout << " CUDA device found \n";
      return 1;
    };
    return -1;
  }
};

// double Displace;
// unsigned short minNumPoints, Crow;
// unsigned short Wsize;
// Qtree qtreeIn = NULL;

// gpu_selector selector;
// cpu_selector selector;
//default_selector selector;
//host_selector selector;

// queue q(selector);
queue q(CUDASelector().select_device());

sycl::device qdevice = q.get_device();
sycl::context qcontext = q.get_context();


uint64_t total_num_nodes = 0;

uint64_t cpu_tree_nodes = 0;

typedef struct QtreeG4_t* QtreeG4;

QtreeG4 array_for_copy = NULL;
int32_t copied_nodes = 0;
int32_t cuda_copied_nodes = 0;

uint64_t Npoints = 0;
// Lpoint** array_all_points = NULL;
Lpoint* array_all_points = NULL;
// static constexpr Lpoint* array_all_points = NULL;

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

  // QtreeG4_t( QtreeG4 p, Vector2D& c, float r, point_alloc_t& alloc ) : 
  //             parent(p), center(c), radius(r), points(alloc) {
                
  //     // quadrants.first = NULL;
  //     // quadrants.second = NULL;
  //     // quadrants.third = NULL;
  //     // quadrants.fourth = NULL;

  //     // quadrants.reserve(4);

  //     // quadrants = static_cast<QtreeG4*>(mallocWrap(4 * sizeof(QtreeG4)));

  //     // for(int i=0; i<4; ++i)
  //     //   quadrants[i] = NULL;

  //     quadrants = NULL;

  //   };
  
  QtreeG4_t( QtreeG4 p, Vector2D& c, float r ) : 
              parent(p), center(c), radius(r) {

      quadrants = NULL;
      numPts = 0;
      points = NULL;
      min = NULL;

      // for(int i=0; i<4; ++i)
      //   quadrants[i] = NULL;

    };
};

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

QtreeG5 cudaArray = NULL;



// FUNCTIONS

void* mallocWrap(size_t size);

// template<typename accNode_t, typename accPoint_t>
// Lpoint gpuSearchIndex( Vector2D& point, accNode_t aNodes, accPoint_t aPoints, float radius, int& numNeighs);


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

// Make a box with center the point and the specified radius
void makeBox(Vector2D *point, float radius, Vector2D *min, Vector2D *max)
{
    // printf("Radio: %.2f\n", radius);
    // printf("Centro: [ %.2lf, %.2lf]\n",point->x,point->y );
    min->x = point->x - radius;
    min->y = point->y - radius;
    // min->z = point->z - radius;
    max->x = point->x + radius;
    max->y = point->y + radius;
    // max->z = point->z + radius;
}

__device__ void makeBox(Vector2D *point, float radius, Vector2D *min, Vector2D *max)
{
    // printf("Radio: %.2f\n", radius);
    // printf("Centro: [ %.2lf, %.2lf]\n",point->x,point->y );
    min->x = point->x - radius;
    min->y = point->y - radius;
    // min->z = point->z - radius;
    max->x = point->x + radius;
    max->y = point->y + radius;
    // max->z = point->z + radius;
}

void makeBox(Vector2D *point, double radiusX, double radiusY, Vector2D *min, Vector2D *max)
{

    min->x = point->x - radiusX;
    min->y = point->y - radiusY;

    max->x = point->x + radiusX;
    max->y = point->y + radiusY;

}

template<typename pT>
int boxInside2D(Vector2D& boxMin, Vector2D& boxMax, pT qt)
{
    if(qt->center.x + qt->radius > boxMax.x ||
       qt->center.y + qt->radius > boxMax.y)
        return 0;

    if(qt->center.x - qt->radius < boxMin.x ||
       qt->center.y - qt->radius < boxMin.y)
        return 0;

    return 1;
}

template<typename pT>
__device__ int boxInside2D(Vector2D& boxMin, Vector2D& boxMax, pT qt)
{
    if(qt->center.x + qt->radius > boxMax.x ||
       qt->center.y + qt->radius > boxMax.y)
        return 0;

    if(qt->center.x - qt->radius < boxMin.x ||
       qt->center.y - qt->radius < boxMin.y)
        return 0;

    return 1;
}

// int boxInside2D(Vector2D boxMin, Vector2D boxMax, QtreeG qt)
// {
//     if(qt->center.x + qt->radius > boxMax.x ||
//        qt->center.y + qt->radius > boxMax.y)
//         return 0;

//     if(qt->center.x - qt->radius < boxMin.x ||
//        qt->center.y - qt->radius < boxMin.y)
//         return 0;

//     return 1;
// }

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

__device__ int insideBox2D(Lpoint *point, Vector2D min, Vector2D max)
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
// int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, Qtree qt)
// {
//     if(qt->center.x + qt->radius < boxMin.x ||
//        qt->center.y + qt->radius < boxMin.y)
//         return 0;

//     if(qt->center.x - qt->radius > boxMax.x ||
//        qt->center.y - qt->radius > boxMax.y)
//         return 0;

//     return 1;
// }

template<typename pT>
int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, pT qt)
{
    if(qt->center.x + qt->radius < boxMin.x ||
       qt->center.y + qt->radius < boxMin.y)
        return 0;

    if(qt->center.x - qt->radius > boxMax.x ||
       qt->center.y - qt->radius > boxMax.y)
        return 0;

    return 1;
}

template<typename pT>
__device__ int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, pT qt)
{
    if(qt->center.x + qt->radius < boxMin.x ||
       qt->center.y + qt->radius < boxMin.y)
        return 0;

    if(qt->center.x - qt->radius > boxMax.x ||
       qt->center.y - qt->radius > boxMax.y)
        return 0;

    return 1;
}

template<typename pT>
int boxTotalOverlap2D(Vector2D boxMin, Vector2D boxMax, pT qt)
{
    if(qt->center.x + qt->radius < boxMax.x ||
       qt->center.y + qt->radius < boxMax.y)
        return 0;

    if(qt->center.x - qt->radius > boxMin.x ||
       qt->center.y - qt->radius > boxMin.y)
        return 0;

    return 1;
}

// int isLeaf(Qtree qt)
// {
//     return qt->quadrants[0] == NULL;
// }

template<typename pT>
int isLeaf(pT qt)
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

// Find the child corresponding a given point
// int quadrantIdx(Lpoint *point, Qtree qtree)
// {
//     int child = 0;

//     if(point->x >= qtree->center.x) child |= 2;
//     if(point->y >= qtree->center.y) child |= 1;
//     // if(point->z >= qtree->center.z) child |= 1;

//     return child;
// }

template<typename pT>
int quadrantIdx(Lpoint *point, pT qtree)
{
    int child = 0;

    if(point->x >= qtree->center.x) child |= 2;
    if(point->y >= qtree->center.y) child |= 1;
    // if(point->z >= qtree->center.z) child |= 1;

    return child;
}


template<typename pT>
Lpoint findValidMin(pT qtree, Vector2D& boxMin, Vector2D& boxMax, int& numInside)
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

template<typename pT>
Lpoint searchNeighborsMin(Vector2D& point, pT qtree, float radius, int& numNeighs)
{
    Vector2D boxMin, boxMax;

    numNeighs = 0;
    makeBox(&point, radius, &boxMin, &boxMax);

    return findValidMin(qtree, boxMin, boxMax, numNeighs);
}


Lpoint searchOverlap(Vector2D* point, Qtree qtree, double radiusX, double radiusY, int& numNeighs)
{
    Vector2D boxMin, boxMax;

    numNeighs = 0;
    makeBox(point, radiusX, radiusY, &boxMin, &boxMax);

    return findValidMin(qtree, boxMin, boxMax, numNeighs);
}

// Lpoint findValidMin(QtreeG4 qtree, Vector2D& boxMin, Vector2D& boxMax, int& numInside)
// {
//     // Lpoint tmp, min = nomin;
//     Lpoint tmp, min = {0,0.0,0.0,99999.0};

//     // if( qtree->numPts > 0 ) { //isEmpty??

//       if( qtree->quadrants == NULL ) //isLeaf??
//       {

//         /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
//         if(boxInside2D(boxMin, boxMax, qtree)) {

//           if(qtree->min->z < min.z)
//             min = *(qtree->min);
//           numInside += qtree->numPts;

//         /*Si está parcialmente solapada, tengo dos opciones*/
//         } else if(boxOverlap2D(boxMin, boxMax, qtree)) {

//           int N = qtree->numPts;
//           Lpoint* pointer = qtree->points;

//           /*Si el mínimo está dentro solo comparo una vez y cuento */
//           if( insideBox2D(qtree->min, boxMin, boxMax) ) {

//             if(qtree->min->z < min.z)
//               min = *(qtree->min);

//             for(int i=0; i<N; ++i) {
//               if(insideBox2D(&pointer[i], boxMin, boxMax))
//               {
//                 numInside++;
//               }
//             }
//           }
//           /*Si el mínimo está fuera, tengo que comprobarlo todo*/
//           else {
//             for(int i=0; i<N; ++i) {
//               if(insideBox2D(&pointer[i], boxMin, boxMax))
//               {
//                 if (pointer[i].z < min.z) {
//                     min = pointer[i];
//                 }
//                 numInside++;
//               }
//             }
//           }

//         }

//       } else {

//           for(int i = 0; i < 4; i++) {

//               QtreeG4 current = &(qtree->quadrants[i]);
//               // Check
//               if(!boxOverlap2D(boxMin, boxMax, current) || current->numPts == 0 ) { //No solapada o vacia
//                   continue;
//               } 
//               /*si la caja está completamente solapada, capturo el mínmo y el numero de puntos, y continuo*/
//               else if( boxInside2D(boxMin, boxMax, current) ) {
//                 if(current->min->z < min.z)
//                   min = *(current->min);
//                 numInside += current->numPts;
//               }
//               /*pero si solo está parcialmente solapada, profundizo*/
//               else {
//                   tmp = findValidMin(current, boxMin, boxMax, numInside);
//                   if (tmp.z < min.z) {
//                       min = tmp;
//                   }
//               }

//           }

//       }

//     // }// isEmpty

//     return min;
// }

// Lpoint searchNeighborsMin(Vector2D& point, QtreeG4 qtree, float radius, int& numNeighs)
// {
//     Vector2D boxMin, boxMax;

//     numNeighs = 0;
//     makeBox(&point, radius, &boxMin, &boxMax);

//     return findValidMin(qtree, boxMin, boxMax, numNeighs);
// }

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

Qtree createQtreeF(Qtree parent, Vector2D center, float radius)
{
    Qtree qt = new Qtree_t;

    qt->center = center;
    qt->radius = radius;
    qt->parent = parent;

    qt->quadrants.reserve(4);
    
    for( int i = 0; i < 4; i++)
      qt->quadrants[i] = NULL;

    cpu_tree_nodes++;

    return qt;
}


void createQuadrantsF(Qtree qt)
{
    Vector2D newCenter;
    float newRadius = qt->radius * 0.5;

    for( int i = 0; i < 4; i++)
    {
        newCenter = qt->center;
        newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);

        qt->quadrants[i] = createQtreeF(qt, newCenter, newRadius);

    }
    // printf("\n");
}


// Insert a point in the qtree creating the appropiate childs
void insertPointF(Lpoint *point, Qtree qtree, float minRadius)
{
    int idx = 0;

    if(isLeaf(qtree))
    {
        // printf("octante hoja nivel %d\n",nivel);
        if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createQuadrantsF(qtree);
          // fillOctants(qtree);
          idx = quadrantIdx(point, qtree);
          insertPointF(point, qtree->quadrants[idx], minRadius);

        } else {
          qtree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
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

void deleteQtree(Qtree qtree)
{
    if(isLeaf(qtree))
    {
      qtree->points.clear();
      qtree->quadrants.clear();
      // delete(qtree);
    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; i++) {
            // Check
            deleteQtree(qtree->quadrants[i]);
            delete(qtree->quadrants[i]);
            // qtree->quadrants.clear();
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }
        qtree->quadrants.clear();
        // delete(qtree);
    }

    return;
}

template<typename pT>
void deleteQtreeGPU(pT qtree)
{
    if(isLeaf(qtree))
    {
      qtree->points.clear();
      free(qtree->quadrants, qcontext);
    } 
    else 
    {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; ++i) {
            // Check
            deleteQtreeGPU(qtree->quadrants[i]);
            free(qtree->quadrants[i], qcontext);
            // delete(qtree->quadrants[0]);
            // qtree->quadrants.clear();
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }
        free(qtree->quadrants, qcontext);
    }

    return;
}

template<typename pT>
void deleteQtreeGPU2(pT qtree) // Este tengo que utilizarlo cuando creo los nodos contiguos
{
    if(isLeaf(qtree))
    {
      qtree->points.clear();
      free(qtree->quadrants, qcontext);
    } 
    else 
    {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; ++i) {
            // Check
            deleteQtreeGPU2(qtree->quadrants[i]);
            // free(qtree->quadrants[i], qcontext);
            // delete(qtree->quadrants[0]);
            // qtree->quadrants.clear();
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }
        free(qtree->quadrants[0], qcontext);
        free(qtree->quadrants, qcontext);
    }

    return;
}

template<typename pT>
void deleteQtreeGPU3(pT qtree) // Este tengo que utilizarlo cuando creo los nodos contiguos
{
    if(isLeaf(qtree))
    {
      qtree->points.clear();
      free(qtree->quadrants, qcontext);
    } 
    else 
    {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; ++i) {
            // Check
            // printf("octante intermedio %d\n", i);
            // fflush(stdout);
            deleteQtreeGPU3(qtree->quadrants[i]);
            // free(qtree->quadrants[i], qcontext);
            // delete(qtree->quadrants[0]);
            // qtree->quadrants.clear();
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }
        // free(qtree->quadrants[0], qcontext);
        free(qtree->quadrants, qcontext);
    }

    return;
}


std::vector<Qtree>::iterator findIterator2(Qtree aux, Qtree current, int& idx) 
{
  std::vector<Qtree>::iterator iterador = std::find(aux->quadrants.begin(), aux->quadrants.end(), current);

  int distance = std::distance(aux->quadrants.begin(), iterador);

  if(distance<3){

    std::advance(iterador, 1);
    idx = distance+1;

  } else { //termino

    iterador = aux->quadrants.begin();
    std::advance(iterador, 3);
    idx = 4;

  }

  // current = *itt;

  return iterador;
}

std::vector<Qtree>::iterator findIterator(Qtree aux, Qtree current, int& idx) 
{
  std::vector<Qtree>::iterator iterador = aux->quadrants.begin();

  if(aux->quadrants[0] == current){
    std::advance(iterador, 1);
    idx=1;
  }
  else if(aux->quadrants[1] == current){
    std::advance(iterador, 2);
    idx=2;
  }
  else if(aux->quadrants[2] == current){
    std::advance(iterador, 3);
    idx=3;
  }
  else {
    //lo pongo apuntando al último
    std::advance(iterador, 3);
    idx = 4; // termino, pero dejo listo current
  }

  return iterador;
}

/*Este se utiliza cuando NO tengo los cuadrantes consecutivos en memora*/
// QtreeG4 findPosition(QtreeG4 parent, QtreeG4 current, int& idx) 
// {

//   if(parent->quadrants[0] == current){
//     idx=1;
//     return parent->quadrants[1];
//   }
//   else if(parent->quadrants[1] == current){
//     idx=2;
//     return parent->quadrants[2];
//   }
//   else if(parent->quadrants[2] == current){
//     idx=3;
//     return parent->quadrants[3];
//   } 
//   else {
//     idx = 4; // termino, pero dejo listo current
//     return parent->quadrants[3];
//   }

// }

/*Este se utiliza cuando tengos los cuadrantes consecutivos en memoria*/
QtreeG4 findPosition(QtreeG4 parent, QtreeG4 current, int& idx) 
{

  if(parent->quadrants == current){
    // std::advance(iterador, 1);
    idx=1;
    return &(parent->quadrants[1]);
  }
  else if(&(parent->quadrants[1]) == current){
    // std::advance(iterador, 2);
    idx=2;
    return &(parent->quadrants[2]);
  }
  else if(&(parent->quadrants[2]) == current){
    // std::advance(iterador, 3);
    idx=3;
    return &(parent->quadrants[3]);
  } 
  else {
    //lo pongo apuntando al último
    // std::advance(iterador, 3);
    idx = 4; // termino, pero dejo listo current
    return &(parent->quadrants[3]);
  }

}

template<typename accNode_t>
__device__ QtreeG5 findPosition(accNode_t& arrayNodes, QtreeG5 parent, QtreeG5 current, int& idx) 
{

  // int32_t quadrants = parent->quadrants;
  QtreeG5 quadrants = &(arrayNodes[parent->quadrants]);

  if(quadrants == current){
    // std::advance(iterador, 1);
    idx=1;
    return &(quadrants[1]);
  }
  else if(&(quadrants[1]) == current){
    // std::advance(iterador, 2);
    idx=2;
    return &(quadrants[2]);
  }
  else if(&(quadrants[2]) == current){
    // std::advance(iterador, 3);
    idx=3;
    return &(quadrants[3]);
  } 
  else {
    //lo pongo apuntando al último
    // std::advance(iterador, 3);
    idx = 4; // termino, pero dejo listo current
    return &(quadrants[3]);
  }

}


// Lpoint cpuSearchNeighborsMin(Vector2D& point, QtreeG4 qtree, float radius, int& numNeighs)
// {
//     int idx = 0;
//     int numInside = 0;
//     Vector2D boxMin, boxMax;

//     // *numNeighs = 0;
//     makeBox(&point, radius, &boxMin, &boxMax);

//     Lpoint min = {0,0.0,0.0,99999.0};

//     QtreeG4 parent = qtree;
//     QtreeG4 current = qtree->quadrants;
//     // QtreeG4 current = qtree->quadrants[0];

//     while(current != qtree) {


//       if(idx > 3) {

//         current = current->parent;

//         if(current != qtree){

//           parent = current->parent;

//           current = findPosition(parent, current, idx);

//         }
        
//       } else {

//         if( current->quadrants == NULL ) { // isLeaf??
//         // if( current->quadrants[0] == NULL ) { // isLeaf??

//           if( current->numPts > 0 ) { //isEmpty??

//             /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
//             if(boxInside2D(boxMin, boxMax, current)) {

//               // int N = current->numPts;
//               // Lpoint* pointer = current->points;
//               // // for(Lpoint* p : current->points) {
//               // for(int i=0; i<N; ++i) {
//               //   if (pointer[i].z < min.z) {
//               //       min = pointer[i];
//               //   }
//               //   // numInside++;
//               // }
//               // numInside += N;

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
//             current = &(parent->quadrants[idx]);
//             // current = parent->quadrants[idx];

//           }

//         } else {
          
//           if(!boxOverlap2D(boxMin, boxMax, current) || current->min == NULL) { //No solapada o vacia

//             idx++;
//             if(idx < 4) {
//               current = &(parent->quadrants[idx]);
//               // current = parent->quadrants[idx];

//             }
//           }
//           /*si la caja está completamente solapada, capturo el mínmo y el numero de puntos, y continuo*/
//           else if( boxInside2D(boxMin, boxMax, current) ) {

//             if(current->min->z < min.z)
//               min = *(current->min);
//             numInside += current->numPts;

//             idx++;
//             if(idx < 4) {
//               current = &(parent->quadrants[idx]);
//               // current = parent->quadrants[idx];

//             }
//           }
//           /*pero si solo está parcialmente solapada, profundizo*/
//           else {

//             idx = 0;
//             parent = current;
//             current = current->quadrants;
//             // current = current->quadrants[0];

//           }
//         }
//       }

//     }

//     numNeighs = numInside;

//     return min;
// }

template<typename pT>
void stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, pT qtreeIn, Vector2D min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    Vector2D cellCenter;

    Lpoint newmin;

    int cellPoints;

    for(int jj = 0 ; jj < Ccol ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < Crow ; ii++ ){

          // printf("CELDA %d,%d <-------------------------------------------------------\n",jj,ii);

            cellCenter.x = initX + ii*Displace;

            // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
            newmin = cpuSearchNeighborsMin(cellCenter, qtreeIn, Wsize/2, cellPoints);
            // newmin = gpuSearchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
            // newmin = gpuSearchIndex(cellCenter, Wsize*0.5, cellPoints);

            if(cellPoints >= minNumPoints ){
              // printf("minimo: %d", newmin.id);
              minIDs.push_back(newmin.id);
              // minIDs[jj*Crow+ii] = newmin.id;
            }
        }

    }

    return;
}

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

                newmin = searchNeighborsMin(cellCenter, qtreeIn, Wsize/2, cellPoints);

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

          newmin = searchNeighborsMin(cellCenter, qtreeIn, Bsize/2, cellPoints);

          if(cellPoints>0){
            minGridIDs.push_back(newmin.id);
          }
        }
      }
    }

    return;
}

template<typename pT>
std::vector<int> stage3reduce(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
                            pT qtreeIn, Qtree grid, Vector2D min)
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

                    // newmin = searchNeighborsMin(cellCenter, qtreeIn, Bsize/2, cellPoints);
                    // newmin = cpuSearchNeighborsMin(cellCenter, qtreeIn, Bsize/2, cellPoints);
                    newmin = gpuSearchNeighborsMin(cellCenter, qtreeIn, Bsize/2, cellPoints);

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

/** Wrapper to handle malloc() nicely */
// void* mallocWrap(size_t size)
// {
//     // void *ptr = malloc_host(size, q);
//     // void *ptr = malloc(size, q, global_alloc);
//     void *ptr = malloc(size, qdevice, qcontext, global_alloc);
//     // void *ptr = aligned_alloc(128, size, q, global_alloc);
//     if(!ptr)
//     {
//         fprintf(stderr, "Not enough memory\n");
//         exit(EXIT_FAILURE);
//     }
//     return ptr;
// }

void* mallocWrap(size_t size)
{
  void *ptr = std::malloc(size);
  // void *ptr = malloc_host(size, q);
  // void *ptr = cl::sycl::malloc(size, qdevice, qcontext, global_alloc);
  if (ptr)
      return ptr;
  else
      throw std::bad_alloc{};
}


template<typename pT>
void createQuadrantsGPU(pT qt)
{
    Vector2D newCenter;
    float newRadius = qt->radius * 0.5;

    for(int i = 0; i < 4; i++)
    {
        newCenter = qt->center;
        newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);
        // newCenter.z += qt->radius * (i&1 ? 0.5f : -0.5f);
        // printf("(%lf,%lf) ", newCenter.x, newCenter.y);

        // qt->quadrants[i] = new Qtree_t(qt, newCenter, newRadius);
        qt->quadrants[i] = createQtreeGPU(qt, newCenter, newRadius);
    }
    // printf("\n");
}


size_t num_pts_copiados = 0;

/*Con este método de copia tengo los cuadrantes alineados*/
template<typename cpunode_t>
void copyQtree2(cpunode_t qt, QtreeG4 parent)
{
  if(isLeaf(qt)) {

    // si este nodo es hoja, solo tengo que copiar los puntos ya que lo he guardado anteriormente
    // quadrants lo dejo a NULL porque es hoja
    // parent->points.insert(parent->points.end(), qt->points.begin(), qt->points.end());
    size_t vsize = qt->points.size();
    if(vsize > 0) {
      parent->points = &(array_all_points[num_pts_copiados]);
      // size_t vsize = qt->points.size();
      parent->numPts = vsize;

      // std::memcpy(&(array_all_points[num_pts_copiados]), qt->points.data(), parent->numPts*sizeof(Lpoint*));
      // num_pts_copiados += parent->numPts;

      Lpoint* min = NULL;
      for(Lpoint* p : qt->points) {
        array_all_points[num_pts_copiados] = *p;
        // num_pts_copiados++;
      
        // Lpoint** min = &(parent->min);
        if(min == NULL) {
          // min = p;
          min = &(array_all_points[num_pts_copiados]);
        } else if(p->z < min->z) {
          min = &(array_all_points[num_pts_copiados]);
        }

        num_pts_copiados++;
      }
      parent->min = min;
    }

  } else {

    // ahora es cuando se que quadrants debe apuntar a algo
    parent->quadrants = &(array_for_copy[copied_nodes]);

    // copio los 4 cuadrantes contiguamente
    for(int i=0; i<4; ++i) {

      cpunode_t src = qt->quadrants[i];

      QtreeG4 dst = &(parent->quadrants[i]);

      // new (dst) QtreeG4_t(  parent, 
      //                       src->center, 
      //                       src->radius, 
      //                       p_alloc );
      new (dst) QtreeG4_t(  parent, 
                            src->center, 
                            src->radius );
      
      // dst->parent = parent;
      // dst->quadrants = NULL;
      // dst->center = src->center;
      // dst->radius = src->radius;
      // dst->numPts = 0;
      // dst->points = NULL;
      // dst->min = NULL;

      copied_nodes++;
    }

    // continúo recursivamente para cada uno de los cuadrantes
    Lpoint* min = NULL;
    for(int i=0; i<4; ++i) {

      QtreeG4 dst = &(parent->quadrants[i]);

      copyQtree2(qt->quadrants[i], dst);

      /* Aquí ya tengo un cuadrante correctamente actualizado, 
      por tanto, cuento los puntos bajo su area y búsco el mínimo */
      // if(dst->quadrants != NULL) {
      //   Lpoint** min = &(parent->min);
      //   for(int j=0; j<4; ++j) {
      //     dst->numPts += dst->quadrants[j].numPts;

      //     if(*min == NULL) {
      //       *min = dst->quadrants[j].min;
      //     } else if(dst->quadrants[j].min->z < (*min)->z) {
      //       *min = dst->quadrants[j].min;
      //     }
      //   }
      // }

      parent->numPts += dst->numPts;

      // if(dst->min == NULL) printf("\t\t Problema\n");

      if(dst->min != NULL) {
        if(min == NULL) {
          min = dst->min;
        } else if(dst->min->z < min->z) {
          min = dst->min;
        }
      }

    }
    parent->min = min;
    // printf("   Termino de analizar en %d\n", nivel);

  }

}

int32_t cuda_num_pts_copiados;

template<typename cpunode_t>
void cudaCopyQtree(cpunode_t qt, int32_t parentID)
{
  QtreeG5 parent = &(cudaArray[parentID]);
  
  if(isLeaf(qt)) {

    // si este nodo es hoja, solo tengo que copiar los puntos ya que lo he guardado anteriormente
    // quadrants lo dejo a NULL porque es hoja
    // parent->points.insert(parent->points.end(), qt->points.begin(), qt->points.end());
    size_t vsize = qt->points.size();
    if(vsize > 0) {
      // parent->points = &(array_all_points[num_pts_copiados]);
      parent->points = cuda_num_pts_copiados;
      // size_t vsize = qt->points.size();
      parent->numPts = vsize;

      // std::memcpy(&(array_all_points[cuda_num_pts_copiados]), qt->points.data(), cudaArray[parent].numPts*sizeof(Lpoint*));
      // cuda_num_pts_copiados += cudaArray[parent].numPts;

      int32_t min = -1;
      for(Lpoint* p : qt->points) {
        array_all_points[cuda_num_pts_copiados] = *p;
        // cuda_num_pts_copiados++;
      
        // Lpoint** min = &(cudaArray[parent].min);
        if(min < 0) {
          // min = p;
          // min = &(array_all_points[cuda_num_pts_copiados]);
          min = cuda_num_pts_copiados;
        // } else if(p->z < min->z) {
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

      QtreeG5 dst = &(cudaArray[quadrantsID + i]);

      new (dst) QtreeG5_t(  parentID, 
                            src->center, 
                            src->radius );
      

      cuda_copied_nodes++;
    }

    // continúo recursivamente para cada uno de los cuadrantes
    int32_t min = -1;
    for(int i=0; i<4; ++i) {

      QtreeG5 dst = &(cudaArray[quadrantsID + i]);

      cudaCopyQtree(qt->quadrants[i], quadrantsID + i);

      parent->numPts += dst->numPts;

      // if(dst->min == NULL) printf("\t\t Problema\n");

      if(dst->min > -1) {
        if(min < 0) {
          min = dst->min;
        } else if(array_all_points[dst->min].z < array_all_points[min].z) {
          min = dst->min;
        }
      }

    }
    parent->min = min;
    // printf("   Termino de analizar en %d\n", nivel);

  }

}

/*Con este método de copia tengo los cuadrantes alineados*/
template<typename cpunode_t>
void cudaQtree(cpunode_t qt)
{
  QtreeG5 root = &(cudaArray[0]);

  new (root) QtreeG5_t( 0, 
                        qt->center, 
                        qt->radius );

  cuda_copied_nodes++;

  // hasta aquí ya tengo copiado root
  cudaCopyQtree(qt, 0);
  
  return;
  
}


/*Con este método de copia tengo los cuadrantes alineados*/
template<typename cpunode_t>
QtreeG4 launchCopy(cpunode_t qt)
{
  QtreeG4 root = &(array_for_copy[0]);

  // new (root) QtreeG4_t(  NULL, 
  //                       qt->center, 
  //                       qt->radius, 
  //                       p_alloc );
  new (root) QtreeG4_t(  NULL, 
                        qt->center, 
                        qt->radius );

  // root->parent = NULL;
  // root->quadrants = NULL;
  // root->center = qt->center;
  // root->radius = qt->radius;
  // root->numPts = 0;
  // root->points = NULL;
  // root->min = NULL;

  copied_nodes++;

  // hasta aquí ya tengo copiado root

  copyQtree2(qt, root);

  // cuento todos los puntos de los cuadrantes y busco el mínimo
  // Lpoint* min = root->min;
  // for(int j=0; j<4; ++j) {
  //   root->numPts += root->quadrants[j].numPts;
    // if(root->min == NULL) {
    //   root->min = root->quadrants[j].min;
    // } else if(root->quadrants[j].min->z < root->min->z) {
    //   root->min = root->quadrants[j].min;;
    // }
  // }
  
  return root;
  
}

/*Con este método de copia lo que hago es alinear los nodos que voy a ir
consultando conforme avanza el método de búsqueda de forma natural*/
// template<typename cpunode_t>
// QtreeG4 copyQtreeP(cpunode_t qt, QtreeG4 parent)
// {
//   QtreeG4 pos;

//   if(isLeaf(qt)) {

//     /*Si este nodo de CPU es hoja, tengo que guardar el nodo y copiar los puntos*/
//     pos = &(array_for_copy[copied_nodes]);

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

//     pos = &(array_for_copy[copied_nodes]);

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



template<typename pT>
pT findPosition(pT parent, pT current, int& idx) 
{

  if(parent->quadrants[0] == current){
    // std::advance(iterador, 1);
    idx=1;
    return parent->quadrants[1];
  }
  else if(parent->quadrants[1] == current){
    // std::advance(iterador, 2);
    idx=2;
    return parent->quadrants[2];
  }
  else if(parent->quadrants[2] == current){
    // std::advance(iterador, 3);
    idx=3;
    return parent->quadrants[3];
  } 
  else {
    //lo pongo apuntando al último
    // std::advance(iterador, 3);
    idx = 4; // termino, pero dejo listo current
    return parent->quadrants[3];
  }

}

/*Este se utiliza cuando NO tengo los cuadrantes consecuivos en memoria*/
// Lpoint gpuSearchNeighborsMin(Vector2D& point, QtreeG4 qtree, float radius, int& numNeighs)
// {
//     int idx = 0;
//     int numInside = 0;
//     Vector2D boxMin, boxMax;

//     // *numNeighs = 0;
//     makeBox(&point, radius, &boxMin, &boxMax);

//     Lpoint min = {0,0.0,0.0,99999.0};

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

/*Este se utiliza cuando tengo todos los cuadrantes consecutivos*/
Lpoint gpuSearchNeighborsMin(Vector2D& point, QtreeG4 qtree, float radius, int& numNeighs)
{
    int idx = 0;
    int numInside = 0;
    Vector2D boxMin, boxMax;

    // *numNeighs = 0;
    makeBox(&point, radius, &boxMin, &boxMax);

    Lpoint min = {0,0.0,0.0,99999.0};

    QtreeG4 parent = qtree;
    QtreeG4 current = qtree->quadrants;
    // QtreeG4 current = qtree->quadrants[0];

    while(current != qtree) {


      if(idx > 3) {

        current = current->parent;

        if(current != qtree){

          parent = current->parent;

          current = findPosition(parent, current, idx);

        }
        
      } else {

        if( current->quadrants == NULL ) { // isLeaf??
        // if( current->quadrants[0] == NULL ) { // isLeaf??

          if( current->numPts > 0 ) { //isEmpty??

            /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
            if(boxInside2D(boxMin, boxMax, current)) {

              // int N = current->numPts;
              // Lpoint* pointer = current->points;
              // // for(Lpoint* p : current->points) {
              // for(int i=0; i<N; ++i) {
              //   if (pointer[i].z < min.z) {
              //       min = pointer[i];
              //   }
              //   // numInside++;
              // }
              // numInside += N;

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

          }

          idx++;
          if(idx < 4) {
            current = &(parent->quadrants[idx]);
            // current = parent->quadrants[idx];

          }

        } else {
          
          if(!boxOverlap2D(boxMin, boxMax, current) || current->min == NULL) { //No solapada o vacia

            idx++;
            if(idx < 4) {
              current = &(parent->quadrants[idx]);
              // current = parent->quadrants[idx];

            }
          }
          /*si la caja está completamente solapada, capturo el mínmo y el numero de puntos, y continuo*/
          else if( boxInside2D(boxMin, boxMax, current) ) {

            if(current->min->z < min.z)
              min = *(current->min);
            numInside += current->numPts;

            idx++;
            if(idx < 4) {
              current = &(parent->quadrants[idx]);
              // current = parent->quadrants[idx];

            }
          }
          /*pero si solo está parcialmente solapada, profundizo*/
          else {

            idx = 0;
            parent = current;
            current = current->quadrants;
            // current = current->quadrants[0];

          }
        }
      }

    }

    numNeighs = numInside;

    return min;
}

template<typename accNode_t, typename accPoint_t>
__device__ Lpoint gpuSearchIndex( Vector2D& point, accNode_t aNodes, accPoint_t aPoints, float radius, int& numNeighs)
{
    int idx = 0;
    int numInside = 0;
    Vector2D boxMin, boxMax;

    // *numNeighs = 0;
    makeBox(&point, radius, &boxMin, &boxMax);

    Lpoint min = {0,0.0,0.0,99999.0};

    QtreeG5 qtree = &(aNodes[0]);

    QtreeG5 parent = qtree;
    QtreeG5 current = &(aNodes[qtree->quadrants]);

    // QtreeG4 parent = qtree;
    // QtreeG4 current = qtree->quadrants;
    // QtreeG4 current = qtree->quadrants[0];

    // printf("Empiezo: \n");

    while(current != qtree) {

      // printf("Empiezo: \n");


      if(idx > 3) {

        // printf(" vuelvo\n");
        
        // current = current->parent;
        current = &(aNodes[current->parent]);

        if(current != qtree){

          // parent = current->parent;
          parent = &(aNodes[current->parent]);

          current = findPosition(aNodes, parent, current, idx);

        }
        
      } else {

        if( current->quadrants < 0 ) { // isLeaf??
        // if( current->quadrants[0] == NULL ) { // isLeaf??
          // printf(" Es hoja\n");

          if( current->numPts > 0 ) { //isEmpty??

            /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
            if(boxInside2D(boxMin, boxMax, current)) {

              // int N = current->numPts;
              // Lpoint* pointer = current->points;
              // // for(Lpoint* p : current->points) {
              // for(int i=0; i<N; ++i) {
              //   if (pointer[i].z < min.z) {
              //       min = pointer[i];
              //   }
              //   // numInside++;
              // }
              // numInside += N;

              Lpoint* currentMin = &(aPoints[current->min]);

              if(currentMin->z < min.z)
                min = *currentMin;
              numInside += current->numPts;

            /*Si está parcialmente solapada, tengo dos opciones*/
            } else if(boxOverlap2D(boxMin, boxMax, current)) {

              int N = current->numPts;
              // Lpoint* pointer = current->points;
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
            // current = &(parent->quadrants[idx]);
            // current = parent->quadrants[idx];
            current = &(aNodes[parent->quadrants + idx]);

          }

        } else {
          
          if(!boxOverlap2D(boxMin, boxMax, current) || current->min < 0) { //No solapada o vacia

            // if(current->min < 0)
              // printf("  %d ", current->min);
              
            // printf("  No solapada o vacia\n");

            idx++;
            if(idx < 4) {
              current = &(aNodes[parent->quadrants + idx]);
              // current = parent->quadrants[idx];

            }
          }
          /*si la caja está completamente solapada, capturo el mínmo y el numero de puntos, y continuo*/
          else if( boxInside2D(boxMin, boxMax, current) ) {

            Lpoint* currentMin = &(aPoints[current->min]);

            if(currentMin->z < min.z)
              min = *currentMin;
            numInside += current->numPts;

            idx++;
            if(idx < 4) {
              current = &(aNodes[parent->quadrants + idx]);
              // current = parent->quadrants[idx];

            }
          }
          /*pero si solo está parcialmente solapada, profundizo*/
          else {
            // printf(" profundizo\n");

            idx = 0;
            parent = current;
            // current = current->quadrants;
            // current = current->quadrants[0];
            current = &(aNodes[current->quadrants]);            

          }
        }
      }

    }

    numNeighs = numInside;

    return min;
}



template<typename pT>
void stage1gpu(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, pT qtree, Vector2D min)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, global_alloc>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    std::cout << " Using device: " << q.get_device().get_info<info::device::name>() << std::endl;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // std::vector<int> values(Ccol);
    // std::iota(std::begin(values), std::end(values), 0);

    q.parallel_for<class minSearch>(range<2>(Ccol,Crow), [=](id<2> index) { 
    // q.parallel_for(range<1>(Ccol), [=](id<1> j) { 

      int jj = static_cast<int>(index[0]);
      int ii = static_cast<int>(index[1]);

      Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

      int cellPoints;

      Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);


      if( cellPoints >= minNumPoints ){
        minGPU[jj*Crow+ii] = newmin.id;
      }

    }).wait();

    // std::vector<int> prueba;

    // prueba.reserve(100);

    // buffer<int, 1> bufout_vect(prueba.data(), range<1>(100));

    // q.submit([&](handler &h) {
    //   // # setup sycl stream class to print standard output from device code
    //   // auto out = stream(4096, 1024, h);

    //   // auto V3 = bufout_vect.get_access<access::mode::write>(h);

    //   // # nd-range kernel
    //   h.parallel_for<class minSearch>(range<2>(Ccol,Crow), [=](id<2> index) { 
    //   // h.parallel_for(range<1>(Ccol), [=](id<1> j) { 

    //     int jj = static_cast<int>(index[0]);
    //     int ii = static_cast<int>(index[1]);

    //     Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

    //     int cellPoints;

    //     Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);

    //     // out << " (" << Ccol << ", " << Crow << ") " <<  "minimo: (" << ii << ", " << jj << ") " << endl;

    //     if( cellPoints >= minNumPoints ){

    //         // out << "minimo: (" << ii << ", " << jj << ") " << endl;
    //         // out << "minimo: " << newmin.id << endl;
    //         minGPU[jj*Crow+ii] = newmin.id;

    //         // v.push_back(newmin.id);
    //         // printf("Tengo un minimo!\n");
    //         // std::cout << "Tengo un minimo\n";

    //     }

    //   });

    // }).wait();


    return;
}

void stage1index(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, Vector2D min)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, global_alloc>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    std::cout << " Using device: " << q.get_device().get_info<info::device::name>() << std::endl;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    auto nR = range<1>(cuda_copied_nodes);
    auto pR = range<1>(cuda_num_pts_copiados);
    auto mR = range<1>(Ccol*Crow);

    buffer<QtreeG5_t, 1> bufNodes(cudaArray, nR);
    buffer<Lpoint, 1> bufPoints(array_all_points, pR);
    buffer<int, 1> bufMins(minGPU, mR);

    q.submit([&](handler &h) {

      /*Hay que poner modo read_write porque de no ser así da un error: 
      " cannot initialize a variable of type 'const QtreeG5' (aka 'QtreeG5_t *const') 
      with an rvalue of type 'const QtreeG5_t *' " cuando se llama a gpuSearchIndex */

      auto accNodes = bufNodes.get_access<access::mode::read_write>(h);
      auto accPoints = bufPoints.get_access<access::mode::read_write>(h);
      auto accC = bufMins.get_access<access::mode::write>(h);

      // sycl::stream out(1024, 256, h);

      // # nd-range kernel
      h.parallel_for<class minSearch>(range<2>(Ccol,Crow), [=](id<2> index) { 
      // h.parallel_for(range<1>(Ccol), [=](id<1> j) { 

        int jj = static_cast<int>(index[0]);
        int ii = static_cast<int>(index[1]);

        Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

        int cellPoints;

        // Lpoint newmin = gpuSearchIndex(cellCenter, accNodes, accPoints, Wsize*0.5, cellPoints);

        // // out << " (" << Ccol << ", " << Crow << ") " <<  "minimo: (" << ii << ", " << jj << ") " << endl;

        // if( cellPoints >= minNumPoints ){

        //     // out << "minimo: (" << ii << ", " << jj << ") " << endl;
        //     // out << "minimo: " << newmin.id << endl;
        //     // minGPU[jj*Crow+ii] = newmin.id;
        //     accC[jj*Crow+ii] = newmin.id;

        //     // v.push_back(newmin.id);
        //     // printf("Tengo un minimo!\n");
        //     // std::cout << "Tengo un minimo\n";

        // }

      });

    // }).wait();
    });
    // bloqueante; ya deberíane estar los datos en minGPU
    auto hMins = bufMins.get_access<access::mode::read>();

    // std::memcpy(minGPU, bufMins, Ccol*Crow*sizeof(int));

    return;
}

/*Este kernel no puede estar descomentado si quiero compilar el cpp*/

// CUDA kernel. Each thread takes care of one element of c
__global__ void vecAdd(int *a, int n) {
  // Get our global thread ID
  int id = blockIdx.x * blockDim.x + threadIdx.x;

  // Make sure we do not go out of bounds
  if (id < n) {
    a[id] = a[id] + 1;
  }
}

template<typename accNode_t, typename accPoint_t, typename accMin_t>
__global__ void cudaSearchKernel(unsigned short Wsize, double Displace, double initX, double initY,
 unsigned short Crow, unsigned short Ccol, unsigned short minNumPoints, accNode_t aNodes, accPoint_t aPoints, accMin_t minIDs) 
 {

    /*En este kernel hay que lanzar 1 solo bloque con 1 solo thread*/
    Vector2D cellCenter;

    Lpoint newmin;

    int cellPoints;

    for(int jj = 0 ; jj < Ccol ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < Crow ; ii++ ){

            cellCenter.x = initX + ii*Displace;

            newmin = gpuSearchIndex(cellCenter, aNodes, aPoints, Wsize*0.5, cellPoints);

            if(cellPoints >= minNumPoints ){
              // printf("minimo %d\n", newmin.id);
              minIDs[jj*Crow+ii] = newmin.id;
            }
        }

    }

    /*Con este kernel tengo un thread por bloque, y analizo un celda por bloque:
      Tengo Ncells bloques de 1 thread*/
    // Vector2D cellCenter = {initX + blockIdx.x * Displace, initY + blockIdx.y * Displace};

    // int cellPoints;

    // Lpoint newmin = gpuSearchIndex(cellCenter, aNodes, aPoints, Wsize*0.5, cellPoints);

    // if(cellPoints >= minNumPoints ){
    //   // printf("minimo %d\n", newmin.id);
    //   minIDs[blockIdx.y * Crow + blockIdx.x] = newmin.id;
    // }

    /*Con este puedo elegir un número arbitrario de bloques e hilos, pero debe haber más de Ncells hilos*/
    // int posX = blockIdx.x * blockDim.x + threadIdx.x;
    // int posY = blockIdx.y * blockDim.y + threadIdx.y;

    // if(posX < Crow && posY < Ccol) {

    //   Vector2D cellCenter = {initX + posX * Displace, initY + posY * Displace};

    //   int cellPoints;

    //   Lpoint newmin = gpuSearchIndex(cellCenter, aNodes, aPoints, Wsize*0.5, cellPoints);

    //   if(cellPoints >= minNumPoints ){
    //     // printf("minimo %d\n", newmin.id);
    //     minIDs[posY * Crow + posX] = newmin.id;
    //   }

    // }

    int countGPU=0;
    for(int i=0; i < Crow*Ccol; ++i){
      // if(minGPU[i]>0) printf("%d\n", minGPU[i]);
      if(minIDs[i]>0){
        countGPU++;
        printf("minimo: %d\n", minIDs[i]);
      }
    }
    printf("\nCeldas no descartadas dentro del kernel CUDA:\t%d\n", countGPU);

}

template<typename pT>
void cudaStage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, pT qtree, Vector2D min)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, global_alloc>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    std::cout << " Using device: " << q.get_device().get_info<info::device::name>() << std::endl;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    auto nR = range<1>(cuda_copied_nodes);
    auto pR = range<1>(cuda_num_pts_copiados);
    auto mR = range<1>(Ccol*Crow);

    buffer<QtreeG5_t, 1> bufNodes(cudaArray, nR);
    buffer<Lpoint, 1> bufPoints(array_all_points, pR);
    buffer<int, 1> bufMins(minGPU, mR);

    q.submit([&](handler &h) {

      auto accA = bufNodes.get_access<access::mode::read_write>(h);
      auto accB = bufPoints.get_access<access::mode::read_write>(h);
      auto accC = bufMins.get_access<access::mode::write>(h);

      // # nd-range kernel
      // h.parallel_for<class minSearch>(range<2>(Ccol,Crow), [Displace,initX,initY](id<2> index) { 
      h.interop_task([=](interop_handler ih) { 
      // h.parallel_for(range<1>(Ccol), [=](id<1> j) { 

        auto dNodes = reinterpret_cast<QtreeG5>(ih.get_mem<backend::cuda>(accA));
        auto dPoints = reinterpret_cast<Lpoint*>(ih.get_mem<backend::cuda>(accB));
        auto dMins = reinterpret_cast<int*>(ih.get_mem<backend::cuda>(accC));

        // vecAdd<<<1, 1>>>(dMins, Crow*Ccol);

        int GridSize = 1;
        // dim3 GridSize = dim3(Crow,Ccol);
        int BlockSize = 1;
        // dim3 BlockSize = dim3(8,8);

        cudaSearchKernel<<< GridSize, BlockSize >>>(Wsize, Displace, initX, initY, Crow, Ccol, minNumPoints, dNodes, dPoints, dMins);

      });

    // }).wait();
    });

    cudaDeviceSynchronize();

    auto hMins = bufMins.get_access<access::mode::read>();

    // int countGPU=0;
    // for(int i=0; i < Crow*Ccol; ++i){
    //   // if(minGPU[i]>0) printf("%d\n", minGPU[i]);
    //   if(minGPU[i]>0){
    //     countGPU++;
    //     printf("minimo: %d\n", minGPU[i]);
    //   }
    // }
    // printf("\nCeldas no descartadas en HOST:\t%d\n", countGPU);

    return;
}

template<typename pT>
void stage1heter(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, pT qtree, Vector2D min, float rate)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, global_alloc>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    // std::cout << " Using device: " << qdevice.get_info<info::device::name>() << std::endl;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // std::vector<int> values(Ccol);
    // std::iota(std::begin(values), std::end(values), 0);

    /*Porcentaje de trabajo para GPU*/
    // double rate = 0.2;
    int chunk = static_cast<int>(Crow * rate);

    // q.prefetch(array_for_copy, cpu_tree_nodes * sizeof(QtreeG4_t));
    // q.prefetch(array_all_points, Npoints * sizeof(Lpoint));

    sycl::event e = q.parallel_for<class minSearch>(range<2>(Ccol, chunk), [=](cl::sycl::id<2> index) { 
    // q.parallel_for(range<1>(Ccol), [=](id<1> j) { 

      int jj = static_cast<int>(index[0]);
      int ii = static_cast<int>(index[1]);

      Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

      int cellPoints;

      Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);


      if( cellPoints >= minNumPoints ){
        minGPU[jj*Crow+ii] = newmin.id;
      }

    });

    tbb::parallel_for( tbb::blocked_range2d<int,int>{0,Ccol,chunk,Crow},
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
                Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);
                // newmin = searchNeighborsMin(cellCenter, cpuqt, Wsize/2, cellPoints);
                // Lpoint newmin = searchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);

                if( cellPoints >= minNumPoints ){

                    minGPU[jj*Crow+ii] = newmin.id;

                }
            }
        }

    });

    e.wait();

    return;
}


template<typename pT>
void stage1acc(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, pT qtree, Vector2D min, float rate)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    /*Porcentaje de trabajo para GPU*/
    // double rate = 0.2;
    int chunk = static_cast<int>(Crow * rate);

    buffer<QtreeG4_t, 1> buffer_nodes(array_for_copy, cpu_tree_nodes);
    buffer<Lpoint, 1> buffer_points(array_all_points, Npoints);
    buffer<int, 1> buffer_mins(minGPU, Crow*Ccol);

    sycl::event e = q.submit([&](handler &h) {

      auto a = buffer_nodes.get_access<access::mode::read>(h);
      auto b = buffer_points.get_access<access::mode::read>(h);
      auto m = buffer_mins.get_access<access::mode::write>(h);
    
      h.parallel_for(range<2>(Ccol, chunk), [=](id<2> index) { 
      // q.parallel_for(range<1>(Ccol), [=](id<1> j) { 

        int jj = static_cast<int>(index[0]);
        int ii = static_cast<int>(index[1]);

        Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

        int cellPoints;

        Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);


        if( cellPoints >= minNumPoints ){
          m[jj*Crow+ii] = newmin.id;
        }

      });

    });

    tbb::parallel_for( tbb::blocked_range2d<int,int>{0,Ccol,chunk,Crow},
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
                Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);
                // newmin = searchNeighborsMin(cellCenter, cpuqt, Wsize/2, cellPoints);
                // Lpoint newmin = searchNeighborsMin(cellCenter, qtree, Wsize*0.5, cellPoints);

                if( cellPoints >= minNumPoints ){

                    minGPU[jj*Crow+ii] = newmin.id;

                }
            }
        }

    });

    e.wait();

    return;
}


unsigned int stage2GPU(unsigned int countMin, int* minIDs){

    unsigned int index = 0;

    int ii,jj,id;

    // for(int i=0; i< countMin; ++i) {
    //   printf("%d ", minIDs[i]);
    // }
    // printf("\n");
    // id = minIDs[0];
    for( jj=0 ;  minIDs[jj] == -1 ; ++jj );

    for( ii=jj ; ii<countMin ; ii=jj ){
        id = minIDs[ii];

        for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

        /* esta es la forma de descartar las posiciones no utilizadas, inicializadas a -1 */
        if(jj-ii > 1){
            // if(jj-ii > 1){
            minIDs[index]=id;
            index++;
            // }
        }
    }

    return index;
}

#endif // ENVI_GPU_HPP