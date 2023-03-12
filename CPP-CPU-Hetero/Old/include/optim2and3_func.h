#ifndef O2AND3_H

#define O2AND3_H

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

#include "tbb/tbb.h"
// #include "tbb/tbbmalloc_proxy.h"
// #include <atomic>

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

using Qtree = struct Qtree_t*;

struct Qtree_t {
  Qtree quadrants[4];
  Vector2D center;
  float radius;
  std::vector<Lpoint*> points;

  Qtree_t(Vector2D c, float r);
};


double round2d(double z);

int check_results(std::string filename, std::vector<int>& ids, Lpoint** pointer, float Displace);

int save_file(std::string file_name, std::vector<int>& ids, Lpoint** pointer);

uint64_t read_points(std::string filename, Lpoint** point_cloud);

int readXYZfile(std::string filename, Lpoint* & point_cloud, unsigned int &Npoints, Vector2D &min, Vector2D &max);

Vector2D getRadius(Vector2D min, Vector2D max, float *maxRadius);

Vector2D getCenter(Vector2D min, Vector2D radius);

void insertPointF(Lpoint *point, Qtree qtree, float minRadius);

void insertPointMaxNumber(Lpoint *point, Qtree qtree, float minRadius, int medSize);

void deleteQtree( Qtree qtree );

Lpoint searchNeighborsMin(Vector2D* point, Qtree qtree, float radius, int& numInside);

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

#endif // ENVI_H
