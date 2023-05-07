#ifndef ENVI_H

#define ENVI_H

#include <iostream>
#include <limits>
#include <cmath>
#include <omp.h>
#include <vector>
#include <functional>

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
    // std::vector<Octree_t> children;
    Octree octants[8];
    Vector3D center;
    Lpoint **points;
    int numPts;
    float radius;
};

typedef struct OctreeI_t* OctreeI;

struct OctreeI_t {
    // std::vector<Octree_t> children;
    OctreeI octants[8];
    Vector3D center;
    Lpoint **points;
    int numPts;
    float radius;
};

#define MIN_RADIUS 1.0 //For the discretization

static const int REALLOC_INCREMENT = 256;

static const Lpoint nomin = {0,0,0,std::numeric_limits<double>::max()};

double round2d(double z);

unsigned int findMin(Lpoint** neighbors, unsigned int cellPoints);

unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage1cpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

int stage1tbb(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

unsigned int stage3cpp(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

unsigned int stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

Octree createOctree(Vector3D center, float radius);

int isLeaf(Octree oct);

int isEmpty(Octree oct);

void insertPoint(Lpoint *point, Octree octree);

void insertPointF(Lpoint *point, Octree octree, float minRadius, int nivel);

void insertPointMinRadius(Lpoint *point, Octree octree, float minRadius);

Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs);

Lpoint searchNeighborsMin(Lpoint* point, Octree octree, float radius, int* numNeighs);

void countNeighbors2D(Lpoint* point, Octree octree, float radius, int* numNeighs);

#endif // ENVI_H
