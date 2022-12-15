#ifndef ENVI_H

#define ENVI_H

#include <iostream>
#include <limits>
#include <cmath>
#include <omp.h>
#include <vector>
#include <unistd.h>
#include <functional>
#include <chrono>
#include <string.h>
#include <tbb/task_scheduler_init.h>

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

// using uOctree = std::unique_ptr<struct Octree_t, free_delete>;

// using uOctree = std::unique_ptr<class Octree_t>;

typedef class Octree_t* Octree;

// struct Octree_t {
//     uOctree octants[8];
//     Vector3D center;
//     float radius;
//     std::vector<Lpoint*> points;
// };

class Octree_t {
public:
  Octree_t(Octree p, Vector3D c, float r) : oparent(p), center(c), radius(r) {
    for(int i = 0; i < 8; i++)
      octants[i] = NULL;
  }
  ~Octree_t() {};

  Octree octants[8];
  Octree oparent;
  Vector3D center;
  float radius;
  std::vector<Lpoint*> points;
};

double round2d(double z);

Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

void insertPointF(Lpoint *point, Octree octree, float minRadius);

void deleteOctree( Octree octree );

unsigned int stage1cpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage1cppParent(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

int stage1tbb(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3cpp(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

// void prueba1( uOctree2& octree);

#endif // ENVI_H
