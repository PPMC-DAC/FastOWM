#ifndef O3_FUNC_H

#define O3_FUNC_H

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
using LpointID = uint32_t;

using Qtree = struct Qtree_t* ;

struct Qtree_t {
  Qtree quadrants[4]; // 4 children
  // Qtree oparent;
  Vector2D center; // center of the 2D box/area covered by this node
  float radius; // radius (halfSide) of the 2D box/area represented by this node
  LpointID min; // ID (index in cloud array) of the minimum point (z-value) covered by this node
  uint32_t numPoints; // number of points hanging from this node
  std::vector<LpointID> points; // in leaf-nodes, stores the LiDAR points in the corresponding area

  Qtree_t(Vector2D c, float r): center{c}, radius{r}, min{0}, numPoints{0} 
  {
    for(int i = 0; i < 4; i++) quadrants[i] = nullptr;
  };
};

using QtreeConc = struct QtreeConc_t* ;

struct QtreeConc_t {
  QtreeConc quadrants[4];
  // QtreeConc oparent;
  Vector2D center;
  float radius;
  tbb::concurrent_vector<LpointID> concurrent_points;

  QtreeConc_t(Vector2D c, float r): center{c}, radius{r} {
    for(int i = 0; i < 4; i++) quadrants[i] = nullptr;
  };
};

double round2d(double z);

Vector2D getRadius(Vector2D &min, Vector2D &max, float& maxRadius);
Vector2D getCenter(Vector2D &min, Vector2D &radius);

template<typename Tree>
void deleteQtree(Tree qtree);

void insertPoint(LpointID point, Lpoint *cloud, Qtree qtree, float minRadius);

//This function can be called with node_delimiter as MinRadius or MaxNumber
template<typename tree_policy> // int for MaxNumber or float for MinRadius
Qtree parallel_qtree( int level, Vector2D center, float radius, Lpoint* cloud, int Npoints, tree_policy policy );

void stage1(ushort Wsize, double Overlap, ushort Crow, ushort Ccol,
  ushort minNumPoints, int* minIDs, Lpoint* cloud, Qtree qtreeIn, Vector2D min);

uint32_t stage2(uint32_t countMin, int* minIDs);

std::vector<int> stage3(ushort Bsize, ushort nCols, ushort nRows,
           Lpoint* cloud, Qtree qtreeIn, Qtree grid, Vector2D min);

#include "../src/optim3_func.cpp"

#endif
