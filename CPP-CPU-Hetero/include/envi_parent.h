#ifndef ENVI_QOLD_H

#define ENVI_QOLD_H

#include <iostream>
#include <limits>
#include <cmath>
#include <omp.h>
#include <vector>
#include <unistd.h>
#include <functional>
#include <chrono>
#include <algorithm>
#include <string.h>

// #define MIN_RADIUS 0.10 //For the discretization

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

typedef struct Octree_t* Octree;

struct Octree_t {
  Octree octants[4];
  Octree parent;
  // Vector3D center;
  Vector2D center;
  float radius;
  std::vector<Lpoint*> points;

  // Octree_t(Vector2D c, float r);
  Octree_t(Octree oct, Vector2D c, float r);
  // ~Octree_t();
};

// typedef struct{
//     unsigned int id;
//     double z;
// } zSearch;

// unsigned int findMin(Lpoint** neighbors, double zzmin, unsigned int cellPoints);

double round2d(double z);

// Lpoint createPoint(unsigned int id, double x, double y, double z);

Vector2D getRadius(Vector2D min, Vector2D max, float *maxRadius);

Vector2D getCenter(Vector2D min, Vector2D radius);

void insertPointF(Lpoint *point, Octree octree, float minRadius);

void deleteOctree(Octree octree);

unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector2D min);

unsigned int stage1cppParent(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector2D min);

unsigned int stage1ParentS(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector2D min);

unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector2D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector2D min);

unsigned int stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector2D min);

#endif // ENVI_QOLD_H
