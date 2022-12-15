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

typedef struct Octree_t* Octree;

// using Octreee = std::unique_ptr<struct Octree_t>;

struct Octree_t {
    // std::vector<Octree_t> children;
    Octree octants[8];
    // std::unique_ptr<Octree_t> octants[8];
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

unsigned int stage1cpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3cpp(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

Octree createOctreeF(Vector3D center, float radius);

int isLeaf(Octree oct);

int isEmpty(Octree oct);

void insertPointF(Lpoint *point, Octree octree, float minRadius, int nivel);

Lpoint searchNeighborsMin(Lpoint* point, Octree octree, float radius, int* numNeighs);

void countNeighbors2D(Lpoint* point, Octree octree, float radius, int* numNeighs);

#endif // ENVI_H
