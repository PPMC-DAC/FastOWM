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
#include <memory>

#define MIN_RADIUS 0.10 //For the discretization

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

// typedef struct Octree_t* Octree;
//
// // struct Octree_t {
// //     uOctree octants[8];
// //     Vector3D center;
// //     float radius;
// //     std::vector<Lpoint*> points;
// // };
//
// struct Octree_t {
//   Octree octants[4];
//   // Octree oparent;
//   // Vector3D center;
//   Vector2D center;
//   float radius;
//   std::vector<Lpoint*> points;
//
//   Octree_t(Vector2D c, float r);
//   // ~Octree_t();
// };

using uOctree = std::unique_ptr<class nNode>;
typedef class nNode* Octree;
typedef std::vector<Lpoint*>& vector_t;

class nNode {
public:
  nNode(Vector2D c, float r) : center(c), radius(r) {
    for(int i = 0; i < 4; i++)
      octants[i] = NULL;
  }
  virtual vector_t getVector() = 0;
  ~nNode() {};

  uOctree octants[4];
  Vector2D center;
  float radius;
};

class nInter : public nNode {
public:
  nInter(Vector2D c, float r) : nNode{c,r} {}
  vector_t getVector() {};
  // ~nInter() {};

};

class nLeaf : public nNode {
public:
  nLeaf(Vector2D c, float r) : nNode{c,r} {}
  vector_t getVector() {
    return vpoints;
  }
  // ~nLeaf() {};
private:
  std::vector<Lpoint*> vpoints;
};


unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

unsigned int stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

// unsigned int findMin(Lpoint** neighbors, double zzmin, unsigned int cellPoints);

double round2d(double z);

unsigned int findMin(Lpoint** neighbors, unsigned int cellPoints);

unsigned int findMin2(Lpoint** neighbors, unsigned int cellPoints);


// Lpoint createPoint(unsigned int id, double x, double y, double z);

Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

// Octree createOctree(Vector3D center, float radius);

int isLeaf(Octree oct);

int isEmpty(Octree oct);

void insertPointF(Lpoint *point, uOctree& octree, float minRadius);

Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs);

void deleteOctree(Octree octree);

#endif // ENVI_QOLD_H
