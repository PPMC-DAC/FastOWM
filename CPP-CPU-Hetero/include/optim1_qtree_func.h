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

typedef struct Qtree_t* Qtree;

// struct Qtree_t {
//     uQtree quadrants[8];
//     Vector3D center;
//     float radius;
//     std::vector<Lpoint*> points;
// };

struct Qtree_t {
  Qtree quadrants[4];
  // Qtree oparent;
  // Vector3D center;
  Vector2D center;
  float radius;
  std::vector<Lpoint*> points;

  Qtree_t(Vector2D c, float r);
  // ~Qtree_t();
};

typedef struct{
    unsigned int id;
    double z;
} zSearch;

unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector3D min);

unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector3D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Qtree qtreeIn, Qtree grid, Vector3D min);

unsigned int stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Qtree qtreeIn, Qtree grid, Vector3D min);

// unsigned int findMin(Lpoint** neighbors, double zzmin, unsigned int cellPoints);

double round2d(double z);

unsigned int findMin(Lpoint** neighbors, unsigned int cellPoints);

unsigned int findMin2(Lpoint** neighbors, unsigned int cellPoints);


// Lpoint createPoint(unsigned int id, double x, double y, double z);

Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

Qtree createQtree(Vector3D center, float radius);

int isLeaf(Qtree oct);

int isEmpty(Qtree oct);

void insertPointF(Lpoint *point, Qtree qtree, float minRadius);

Lpoint** searchNeighbors2D(Lpoint *point, Qtree qtree, float radius, int *numNeighs);

void deleteQtree(Qtree qtree);

#endif // ENVI_QOLD_H
