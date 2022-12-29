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
// #include <tbb/task_scheduler_init.h>
#include <memory>

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
// using uOctree = std::unique_ptr<struct Octree_t>;
//
// typedef struct Octree_t* Octree;
//
// struct Octree_t {
//     // std::vector<Octree_t> children;
//     // Octree octants[8];
//     uOctree octants[8];
//     Vector3D center;
//     // Lpoint **points;
//     std::vector<Lpoint*> vpoints;
//     // int numPts;
//     float radius;
// };



// class Node_t {
// public:
//   // virtual bool isInter() const = 0;
//   virtual float myRadius() = 0;
//   float radius;
//   virtual ~Node_t() {};
// };

// using hOctree = std::unique_ptr<Node_t>;

// using uOctree = std::unique_ptr<class nIntermedio>;
// typedef class nIntermedio* Octree;

// class nIntermedio {
// public:
//   // nIntermedio() {};
//   // nIntermedio() {
//     // radius = 0.0;
//     // center = {0.0,0.0,0.0};
//     // for(int i = 0; i < 8; i++)
//     //   octants[i] = NULL;
//   // }
//   nIntermedio(Vector3D c, float r) : radius(r), center(c) {
//     // radius = r;
//     // center = c;
//     for(int i = 0; i < 8; i++)
//       octants[i] = NULL;
//   }
//   // float myRadius() { return radius; }
//   // bool isInter() const override {
//   //   return octants[0] == NULL;
//   // }
//   ~nIntermedio() {};
// // private:
//   uOctree octants[8];
//   Vector3D center;
//   float radius;
//   std::vector<Lpoint*> vpoints;
// };
//
// class nHoja : public nIntermedio {
// public:
//   // nHoja() {};
//   // nHoja() {
//     // radius = 0.0;
//     // center = {0.0,0.0,0.0};
//     // for(int i = 0; i < 8; i++)
//     //   octants[i] = NULL;
//   // }
//   nHoja(Vector3D c, float r) : nIntermedio{c,r} {
//     // radius = r;
//     // center = c;
//     // for(int i = 0; i < 8; i++)
//     //   octants[i] = NULL;
//   }
//   // bool isInter() const override {
//   //   return octants[0] == NULL;
//   // }
//   ~nHoja() {};
// // private:
//   // hOctree octants[8];
//   // std::vector<Lpoint*> vpoints;
//   // Vector3D center;
//   // float radius;
// };



using uOctree = std::unique_ptr<class nNode2>;
typedef class nNode2* Octree;
typedef std::vector<Lpoint*>* vector_t;

class nNode2 {
public:
  nNode2(Vector3D c, float r) : center(c), radius(r) {
    for(int i = 0; i < 8; i++)
      octants[i] = NULL;
  }
  virtual vector_t getVector() = 0;
  ~nNode2() {};

  uOctree octants[8];
  Vector3D center;
  float radius;
};

class nIntermedio2 : public nNode2 {
public:

  nIntermedio2(Vector3D c, float r) : nNode2{c,r} {}
  vector_t getVector() override {};
  ~nIntermedio2() {};

};

class nHoja2 : public nNode2 {
public:
  nHoja2(Vector3D c, float r) : nNode2{c,r} {}
  vector_t getVector() override {
    return &vpoints;
  }
  ~nHoja2() {};
private:
  std::vector<Lpoint*> vpoints;
};



#define MIN_RADIUS 1.0 //For the discretization

static const int REALLOC_INCREMENT = 256;

double round2d(double z);

Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius);

Vector3D getCenter(Vector3D min, Vector3D radius);

// uOctree createOctreeF(Vector3D center, float radius);

// int isLeaf(Octree oct);
// int isLeaf(uOctree& oct);

// int isEmpty(Octree oct);

void insertPointF(Lpoint *point, uOctree& octree, float minRadius, int nivel);
// void insertPointF(Lpoint *point, uOctree* octree, float minRadius, int nivel);

// Lpoint searchNeighborsMin(Lpoint* point, uOctree& octree, float radius, int* numNeighs);
//
// void countNeighbors2D(Lpoint* point, uOctree& octree, float radius, int* numNeighs);

void deleteOctree( Octree octree );

unsigned int stage1cpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min);

unsigned int stage2(unsigned int countMin, int* minIDs);

unsigned int stage3cpp(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min);

// void prueba1( uOctree2& octree);

#endif // ENVI_H
