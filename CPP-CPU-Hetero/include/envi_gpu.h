#ifndef ENVI_H

#define ENVI_H

#pragma once

#include <iostream>
#include <limits>
#include <cmath>
#include <omp.h>
#include <vector>
#include <unistd.h>
#include <functional>
#include <algorithm>
#include <string.h>
#include <numeric>
#include <CL/sycl.hpp>

#include "tbb/tbb.h"
// #include "tbb/tbbmalloc_proxy.h"
// #include <atomic>

// extern int NUM_PROCS;

using namespace cl::sycl;

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

// struct free_delete
// {
//     void operator()(void* x) { free(x); }
// };

// using uQtree = std::unique_ptr<struct Qtree_t, free_delete>;

// using uQtree = std::unique_ptr<class Qtree_t>;

typedef struct Qtree_t* Qtree;

// struct Qtree_t {
//     uQtree quadrants[8];
//     Vector3D center;
//     float radius;
//     std::vector<Lpoint*> points;
// };

struct Qtree_t {
  // Qtree quadrants[4];
  std::vector<Qtree> quadrants;
  Qtree parent;
  // Vector3D center;
  Vector2D center;
  float radius;
  std::vector<Lpoint*> points;

  // Qtree_t(Vector2D c, float r, Qtree parent);
  // ~Qtree_t();
};

typedef struct QtreeG_t* QtreeG;

extern queue q;

using point_alloc_t = usm_allocator<Lpoint*, usm::alloc::host>;

using qtree_alloc_t = usm_allocator<QtreeG_t, usm::alloc::host>;

using addr_alloc_t = usm_allocator<QtreeG, usm::alloc::host>;

using point_vector_t = std::vector<Lpoint*, point_alloc_t>;

using qtree_vector_t = std::vector<QtreeG_t, qtree_alloc_t>;

using addr_vector_t = std::vector<QtreeG, addr_alloc_t>;

extern point_alloc_t p_alloc;

extern qtree_vector_t all_nodes;

extern addr_vector_t all_addresses;

void* mallocWrap(size_t size);

struct QtreeG_t {
  QtreeG* quadrants;
  QtreeG parent;
  // Vector3D center;
  Vector2D center;
  float radius;
  // std::vector<Lpoint*> points;
  point_vector_t points;

  // QtreeG_t(){};
  QtreeG_t( point_alloc_t& alloc ) : points(alloc) {};

  QtreeG_t( QtreeG p, Vector2D c, float r, point_alloc_t& alloc ) : 
              parent(p), center(c), radius(r), points(alloc) {

      quadrants = static_cast<QtreeG*>(mallocWrap(4 * sizeof(QtreeG)));

      for(int i = 0; i < 4; i++)
        quadrants[i] = NULL;
    };
};


double round2d(double z);

int check_results(std::string filename, std::vector<int>& ids, Lpoint** pointer, float Displace);

int save_file(std::string file_name, std::vector<int>& ids, Lpoint** pointer);

uint64_t read_points(std::string filename, Lpoint** point_cloud);

Vector2D getRadius(Vector2D min, Vector2D max, float *maxRadius);

Vector2D getCenter(Vector2D min, Vector2D radius);

void insertPointF(Lpoint *point, Qtree qtree, float minRadius);

// template<typename nT>
// void insertPointGPU(Lpoint *point, nT qtree, float minRadius);

int quadrantIdx(Lpoint *point, QtreeG qtree);

void createQuadrantsGPU3(QtreeG qt);

void insertPointF2(Lpoint *point, Qtree qtree, float minRadius, int medSize);

void deleteQtree( Qtree qtree );

void deleteQtreeGPU(QtreeG qtree);

void deleteQtreeGPU2(QtreeG qtree);
void deleteQtreeGPU3(QtreeG qtree);

Lpoint searchNeighborsMin(Vector2D* point, Qtree qtree, float radius, int* numNeighs);

void stage1gpu(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, QtreeG qtree, Vector2D min);

void stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, Qtree qtreeIn, Vector2D min);

void stage1cpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, Qtree qtreeIn, Vector2D min);

void stage1rem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, Qtree qtreeIn, Vector2D min);

void stage1remCpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, Qtree qtreeIn, Vector2D min);

// unsigned int stage1cppParent(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
//   unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector2D min);

std::vector<int> stage1tbb(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);
std::vector<int> stage1tbbRem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);
std::vector<int> stage1tbb2(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);
std::vector<int> stage1tg(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);
std::vector<int> stage1tgRem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);
std::vector<int> stage1tg2(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);
std::vector<int> stage1task(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);

std::vector<int> stage1reduce(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);
std::vector<int> stage1reduceRem(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
                        unsigned short minNumPoints, Qtree qtreeIn, Vector2D min);

unsigned int stage2(unsigned int countMin, std::vector<int>& minIDs);
unsigned int stage2GPU(unsigned int countMin, int* minIDs);

std::vector<int> stage2cpp(unsigned int countMin, std::vector<int>& minIDs);

void stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          std::vector<int>& minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min);


void stage3cpp(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          std::vector<int>& minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min);

uint32_t stage3tbb(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          std::vector<int>& minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min);

std::vector<int> stage3reduce(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
                            Qtree qtreeIn, Qtree grid, Vector2D min);
// void prueba1( uQtree2& qtree);

Lpoint gpuSearchNeighborsMin(Vector2D* point, QtreeG qtree, float radius, int* numNeighs);

Qtree createQtreeF(Qtree parent, Vector2D center, float radius);

QtreeG createQtreeGPU(QtreeG parent, Vector2D center, float radius);
QtreeG createQtreeGPU2(QtreeG parent, Vector2D center, float radius);

#endif // ENVI_H
