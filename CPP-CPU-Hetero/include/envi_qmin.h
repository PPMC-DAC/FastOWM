#ifndef ENVI_QMIN_H

#define ENVI_QMIN_H

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
    unsigned int id;
    double x;
    double y;
    double z;
} Lpoint;

using Qtree = struct Qtree_t* ;

struct Qtree_t {
  Qtree quadrants[4];
  Vector2D center;
  float radius;
  std::vector<Lpoint*> points;

  Qtree_t(Vector2D c, float r);
};

typedef struct{
    unsigned int id;
    double z;
} zSearch;

double round2d(double z);

Vector2D getRadius(Vector2D &min, Vector2D &max, float *maxRadius);

Vector2D getCenter(Vector2D &min, Vector2D &radius);

void insertPointF(Lpoint *point, Qtree qtree, float minRadius);

void deleteQtree(Qtree qtree);

void stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector2D min);

// unsigned int stage1cppParent(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
//   unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector2D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min);

#endif // ENVI_QMIN_H
