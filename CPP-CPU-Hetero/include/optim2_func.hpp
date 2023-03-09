#ifndef O2_FUNC_H

#define O2_FUNC_H

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>

#include <tbb/tbb.h>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <functional>
#include <chrono>
#include <algorithm>
#include <string.h>
#include <fstream>
#include <sstream>
#include <filesystem>

// Processes input arguments
#include "../include/cxxopts.hpp"

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

using Lpoint = Vector3D;
using LpointID = int;

using Qtree = struct Qtree_t* ;

struct Qtree_t {
  Qtree quadrants[4];
  // Qtree oparent;
  Vector2D center;
  float radius;
  std::vector<LpointID> points;

  Qtree_t(Vector2D c, float r);
};

using QtreeConc = struct QtreeConc_t* ;

struct QtreeConc_t {
  QtreeConc quadrants[4];
  // QtreeConc oparent;
  Vector2D center;
  float radius;
  tbb::concurrent_vector<LpointID> concurrent_points;

  QtreeConc_t(Vector2D c, float r);
};

double round2d(double z);

Vector2D getRadius(Vector2D &min, Vector2D &max, float *maxRadius);

Vector2D getCenter(Vector2D &min, Vector2D &radius);

template<typename Tree>
void deleteQtree(Tree qtree);

void insertPoint(LpointID point, Lpoint *cloud, Qtree qtree, float minRadius);

//This function can be called with node_delimiter as MinRadius or MaxNumber
template<typename tree_policy> // int for MaxNumber or float for MinRadius
Qtree parallel_qtree( int level, Vector2D center, float radius, Lpoint* cloud, int Npoints, tree_policy policy );

void stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Lpoint* cloud, Qtree qtreeIn, Vector2D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

std::vector<int> stage3(unsigned short Bsize, unsigned short nCols, unsigned short nRows,
           Lpoint* cloud, Qtree qtreeIn, Qtree grid, Vector2D min);

#include "../src/optim2_func.cpp"

#endif
