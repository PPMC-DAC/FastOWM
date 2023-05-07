#include "lbvh.cuh"
#include <nppdefs.h>
#include <random>
#include <vector>
#include <chrono>
#include <thrust/random.h>
#include <thrust/pair.h>
// #include <tbb/parallel_invoke.h>

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include "tbb/global_control.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range2d.h"
#include "tbb/concurrent_vector.h"

// #include <cub/cub.cuh>
// #include <cub/device/device_radix_sort.cuh>
#include <bitset>
#include <fstream>


#pragma once

using real_t = double;
using real_v = double2;
using point_t = double4;
// using node_t = lbvh::detail::node;

struct node
{
    std::uint32_t parent_idx; // parent node
    std::uint32_t left_idx;   // index of left  child node
    std::uint32_t right_idx;  // index of right child node
    std::uint32_t numPts;     // number of points
    std::uint32_t min_idx;    // index of subtree's min point
    std::uint32_t object_idx; // == 0xFFFFFFFF if internal node.

};

using node_t = node;

node_t default_node ={ 0xFFFFFFFFu,   // parent
                             0xFFFFFFFFu,   // left
                             0xFFFFFFFFu,   // right
                             0xFFFFFFFFu,   // numPts
                             0xFFFFFFFFu,   // min
                             0xFFFFFFFFu }; // object

struct aabb_t
{
    real_v upper;
    real_v lower;
};

const auto inf = std::numeric_limits<real_t>::infinity();

#define TOL 1e-10
#define FW_S32_MIN  (~0x7FFFFFFF)


template<typename pointer_t>
int isLeaf(pointer_t bt);

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
    uint32_t id;
    double x;
    double y;
    double z;

} Lpoint;


struct tree_t;

typedef struct tree_t* tree;

struct tree_t {
  tree parent;
  tree childs;
  aabb_t bbox;
  int numPts;
  point_t* points;
  point_t* min;
  
  tree_t( tree p, const aabb_t& bb ) : parent(p), bbox(bb) {

      childs = NULL;
      numPts = 0;
      points = NULL;
      min = NULL;

    };
};

struct tree2_t; // GPU node

typedef struct tree2_t* tree2;

struct tree2_t {
  tree2 parent;
  tree2 childs;
  // aabb_t bbox;
  Vector2D center;
  float radius;
  int numPts;
  point_t* points;
  point_t* min;
  
  tree2_t( tree2 p, Vector2D c, float r ) : parent(p), center(c), radius(r) {

      childs = NULL;
      numPts = 0;
      points = NULL;
      min = NULL;

    };
};

struct Btree_t; // GPU node

typedef struct Btree_t* Btree;

struct Btree_t {
  Btree parent;
  Btree childs[2];
  // aabb_t bbox;
  Vector2D center;
  float radius;
  // int numPts;
  std::vector<point_t*> points;
  // point_t* min;
  
  Btree_t( Btree p, Vector2D c, float r ) : parent(p), center(c), radius(r) {

      childs[0] = NULL;
      childs[1] = NULL;
      // numPts = 0;
      // points = NULL;
      // min = NULL;

    };
};

struct LBVH
{
    const node_t*     node_list;
    const aabb_t*     aabb_list;
    const uint32_t*   objectOrder;
    const point_t*    point_list;
    uint32_t          numNodes; /*internal nodes*/

    __device__ __forceinline__ uint32_t  getNumLeaves    (void)           { return numNodes + 1; }
    __device__ __forceinline__ uint32_t  getLeaf         (uint32_t idx)   { return numNodes + idx; }
    __device__ __forceinline__ uint32_t  getRoot         (void)           { return 0; }


    __device__ __forceinline__ bool      isLeaf          (uint32_t node) { return (node >= numNodes); }
    __device__ __forceinline__ uint32_t  getNumPoints    (uint32_t node) { return node_list[node].numPts; }      
    // __device__ __forceinline__ uint32_t  getObjectIdx    (uint32_t node) { return objectOrder[node_list[node].object_idx]; }
    __device__ __forceinline__ const uint32_t* getObjectList   (uint32_t node) { return &(objectOrder[node_list[node].object_idx]); }
    __device__ __forceinline__ uint32_t  getPointIdx    (uint32_t node, uint32_t i) { return objectOrder[node_list[node].object_idx + i]; }
    __device__ __forceinline__ const aabb_t&   getAABB         (uint32_t node) { return aabb_list[node]; }
    __device__ __forceinline__ uint32_t  getLeftChild    (uint32_t node) { return node_list[node].left_idx; }
    __device__ __forceinline__ uint32_t  getRightChild   (uint32_t node) { return node_list[node].right_idx; }
    __device__ __forceinline__ point_t   getMin          (uint32_t node) { return point_list[node_list[node].min_idx]; }
    __device__ __forceinline__ uint32_t  getMinIdx       (uint32_t node) { return node_list[node].min_idx; }
    // __device__ __forceinline__ uint32_t  getRightmostLeafInLeftSubtree  (uint32_t node) { return nodeChildren[node].z + numNodes; }
    // __device__ __forceinline__ uint32_t  getRightmostLeafInRightSubtree (uint32_t node) { return nodeChildren[node].w + numNodes; }
    __device__ __forceinline__ 
    point_t getPoint(uint32_t node, uint32_t i) { 
      point_t aux = point_list[objectOrder[node_list[node].object_idx + i]];
      /*aprovecho el cuarto elemento para enviar el índice*/
      aux.w = real_t(objectOrder[node_list[node].object_idx + i]);
      return aux; 
    }
    __device__ __forceinline__
    thrust::pair<real_t, uint32_t> getPointTuple (uint32_t node, uint32_t i) {
      return thrust::make_pair( point_list[objectOrder[node_list[node].object_idx + i]].z, objectOrder[node_list[node].object_idx + i] );
    }
    __device__ __forceinline__
    thrust::pair<real_t, uint32_t> getMinTuple (uint32_t node) {
      return thrust::make_pair( point_list[node_list[node].min_idx].z, node_list[node].min_idx );
    }
};


/*Esta estructura se utiliza con puntos ordenados*/
struct LBVHo
{
    const node_t*     node_list;
    const aabb_t*           aabb_list;
    const point_t*          point_list;
    // uint32_t                numInternalNodes;

    __device__ __forceinline__ uint32_t  getRoot         (void)           { return 0; }

    /*Ahora hay nodos hoja que no van a cumplir la desigualdad con numInternalNodes*/
    __device__ __forceinline__ bool      isLeaf          (uint32_t node) { return (node_list[node].object_idx != 0xFFFFFFFFu); }
    __device__ __forceinline__ uint32_t  getNumPoints    (uint32_t node) { return node_list[node].numPts; }      
    __device__ __forceinline__ const aabb_t&  getAABB    (uint32_t node) { return aabb_list[node]; }
    __device__ __forceinline__ uint32_t  getLeftChild    (uint32_t node) { return node_list[node].left_idx; }
    __device__ __forceinline__ uint32_t  getRightChild   (uint32_t node) { return node_list[node].right_idx; }

    __device__ __forceinline__ 
    point_t getPoint(uint32_t node, uint32_t i) { 
      int id = node_list[node].object_idx + i;
      point_t aux = point_list[id];
      /*aprovecho el cuarto elemento para enviar el índice*/
      aux.w = real_t(id);
      return aux; 
    }

    __device__ __forceinline__
    thrust::pair<real_t, uint32_t> getMinTuple (uint32_t node) {
        return thrust::make_pair( point_list[node_list[node].min_idx].z, node_list[node].min_idx );
    }
};


extern tree array_nodes;
extern tree2 array_nodes2;
extern point_t* array_all_points;
extern int copied_nodes;
extern int num_pts_copiados;

extern size_t num_objects;
extern uint32_t r_num_objects;
extern uint32_t  r_num_internal_nodes;
extern uint32_t r_num_nodes;



double round2d(double z){
  return round(z*100.0)/100.0;
}

__device__ __host__ __forceinline__
int boxOverlap2D(const aabb_t& global, const aabb_t& local)
{
  // std::cout << global.upper.x << "," << global.upper.y << " ";
  // std::cout << global.lower.x << "," << global.lower.y << "\n";
  // std::cout << local.upper.x << "," << local.upper.y << " ";
  // std::cout << local.lower.x << "," << local.lower.y << "\n\n";

  if((local.upper.x - global.lower.x < TOL) ||
    (local.upper.y - global.lower.y < TOL))
      return 0;

  if((local.lower.x - global.upper.x > TOL) ||
    (local.lower.y - global.upper.y > TOL))
      return 0;

  return 1;

  // return (a.lower.x <= b.upper.x & a.upper.x >= b.lower.x & a.lower.y <= b.upper.y & a.upper.y >= b.lower.y );
}

/* Checks if the search area is completely overlapped by the node */
__device__ __host__ __forceinline__
int isAllOverlaped(const aabb_t& global, const aabb_t& local)
{
  if((local.upper.x - global.upper.x > TOL) ||
    (local.upper.y - global.upper.y > TOL))
      return 0;

  if((local.lower.x - global.lower.x < TOL) ||
    (local.lower.y - global.lower.y < TOL))
      return 0;

  return 1;
}

/* Checks if a point is contained in the search area */
__device__ __host__ __forceinline__
int insideBox(aabb_t& global, point_t& point)
{
    if((point.x - global.lower.x > TOL) && (point.y - global.lower.y > TOL))
    {
        if((point.x - global.upper.x < TOL) && (point.y - global.upper.y < TOL))
        {
            return 1;
        }
    }

    return 0;
}

/* Recursive search of the minimum in an area */
void countPoints(node_t* n, aabb_t* bb, point_t* p, uint32_t* idxs, 
  node_t& current, aabb_t& searchBB, uint32_t& numPts)
{

  if(current.left_idx == 0xFFFFFFFFu)
  {
    const uint32_t* indices = &idxs[current.object_idx];
    uint32_t countPoints=0;
    for(int i=0; i<current.numPts; i++){
      if(insideBox(searchBB, p[indices[i]])){
        countPoints++;
      }
    }
    numPts += countPoints;

  } else {

    int childIdx = current.left_idx;
    if(boxOverlap2D(searchBB, bb[childIdx])){
      if(isAllOverlaped(searchBB, bb[childIdx])){
        numPts += n[childIdx].numPts;
      }
      else {
        countPoints(n,bb,p,idxs,n[childIdx], searchBB, numPts);
      }
    }

    childIdx = current.right_idx;
    if(boxOverlap2D(searchBB, bb[childIdx])){
      if(isAllOverlaped(searchBB, bb[childIdx])){
        numPts += n[childIdx].numPts;
      }
      else {
        countPoints(n,bb,p,idxs,n[childIdx], searchBB, numPts);
      }
    }

  }

  return;
}

// __device__
// point_t traverseRecursive(LBVH& lbvh, aabb_t& queryAABB, uint32_t node)
// {
//     // Bounding box overlaps the query => process node.

//     point_t lmin = {0.0,0.0,9999.0,0.0};
//     uint32_t numPts = 0;

//     if (boxOverlap2D(queryAABB, lbvh.getAABB(node)))
//     {
//         // Leaf node => report collision.

//         if (lbvh.isLeaf(node)){
//             // list.add(queryObjectIdx, lbvh.getObjectIdx(node));
//             const uint32_t* idx = lbvh.getObjectList(node);
//             for(int i=0; i<numPts; i++){
//               point_t p = lbvh.getPoint(node, i);
//               if(insideBox(queryAABB, p)){
//                 numPts++;
//                 if(p.z < lmin.z)
//                   lmin = p;
//               }
//             }
//         }

//         // Internal node => recurse to children.

//         else
//         {
//             // auto mytuple = std::make_tuple(lmin, numPts);
//             uint32_t childL = lbvh.getLeftChild(node);
//             uint32_t childR = lbvh.getRightChild(node);
//             point_t tmp1 = traverseRecursive(lbvh, queryAABB, childL);
//             point_t tmp2 = traverseRecursive(lbvh, queryAABB, childR);
//             if(tmp1.z < lmin.z)
//               lmin = tmp1;
//             if(tmp2.z < lmin.z)
//               lmin = tmp2;
//         }
//     }

//     return lmin;
// }

__device__
uint32_t traverseRecursive(LBVH& lbvh, aabb_t& queryAABB, uint32_t node)
{
    // Bounding box overlaps the query => process node.

    point_t lmin = {0.0,0.0,9999.0,0.0};
    uint32_t numPts = 0;

    if (boxOverlap2D(queryAABB, lbvh.getAABB(node)))
    {
        // Leaf node => report collision.
        if (lbvh.isLeaf(node)){
            // list.add(queryObjectIdx, lbvh.getObjectIdx(node));
            const uint32_t l = lbvh.getNumPoints(node);
            for(int i=0; i<l; i++){
              point_t p = lbvh.getPoint(node, i);
              if(insideBox(queryAABB, p)){
                numPts++;
                if(p.z < lmin.z)
                  lmin = p;
              }
            }
        }

        // Internal node => recurse to children.

        else
        {
            // std::tuple<point_t,uint32_t> mytuple = std::make_tuple(lmin, numPts);
            uint32_t childL = lbvh.getLeftChild(node);
            uint32_t childR = lbvh.getRightChild(node);
            uint32_t tmp1 = traverseRecursive(lbvh, queryAABB, childL);
            uint32_t tmp2 = traverseRecursive(lbvh, queryAABB, childR);
            numPts += tmp1;
            numPts += tmp2;
            // if(tmp1.z < lmin.z)
            //   lmin = tmp1;
            // if(tmp2.z < lmin.z)
            //   lmin = tmp2;
        }
    }

    return numPts;
}

__device__ 
void traverseIterative(LBVH& lbvh, aabb_t& queryAABB, uint32_t& numPts, point_t& lmin)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    uint32_t stack[64];
    uint32_t* stackPtr = stack;
    *stackPtr++ = NULL; // push

    // point_t lmin = {0.0,0.0,9999.0,0.0};
    // uint32_t numPts = 0;


    // Traverse nodes starting from the root.

    uint32_t node = lbvh.getRoot();
    do
    {
        // Check each child node for overlap.

        uint32_t childL = lbvh.getLeftChild(node);
        uint32_t childR = lbvh.getRightChild(node);
        bool overlapL = (boxOverlap2D(queryAABB, lbvh.getAABB(childL)));
        bool overlapR = (boxOverlap2D(queryAABB, lbvh.getAABB(childR)));

        // Query overlaps a leaf node => report collision.

        if (overlapL && lbvh.isLeaf(childL)){
          const uint32_t l = lbvh.getNumPoints(childL);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childL, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < lmin.z)
                lmin = p;
            }
          }
        }
        if (overlapR && lbvh.isLeaf(childR)){
          const uint32_t l = lbvh.getNumPoints(childR);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childR, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < lmin.z)
                lmin = p;
            }
          }
        }

        // Query overlaps an internal node => traverse.

        bool traverseL = (overlapL && !lbvh.isLeaf(childL));
        bool traverseR = (overlapR && !lbvh.isLeaf(childR));

        if (!traverseL && !traverseR)
            node = *--stackPtr; // pop
        else
        {
            node = (traverseL) ? childL : childR;
            if (traverseL && traverseR)
                *stackPtr++ = childR; // push
        }
    }
    while (node != NULL);

    // count = numPts;

    return;
}

__device__ 
void traverseIterative2(LBVH& lbvh, aabb_t& queryAABB, uint32_t& numPts, point_t& lmin)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    uint32_t stack[64];
    uint32_t* stackPtr = stack;
    *stackPtr++ = NULL; // push

    // Traverse nodes starting from the root.

    uint32_t node = lbvh.getRoot();
    do
    {
        // Check each child node for overlap.

        uint32_t childL = lbvh.getLeftChild(node);
        uint32_t childR = lbvh.getRightChild(node);
        bool overlapL = (boxOverlap2D(queryAABB, lbvh.getAABB(childL)));
        bool overlapR = (boxOverlap2D(queryAABB, lbvh.getAABB(childR)));

        // Internal node totally overlaped

        if(isAllOverlaped(queryAABB, lbvh.getAABB(childL))){
          overlapL = false;
          numPts += lbvh.getNumPoints(childL);
          point_t min = lbvh.getMin(childL);
          if(min.z < lmin.z)
            lmin = min;
        }
        if(isAllOverlaped(queryAABB, lbvh.getAABB(childR))){
          overlapR = false;
          numPts += lbvh.getNumPoints(childR);
          point_t min = lbvh.getMin(childR);
          if(min.z < lmin.z)
            lmin = min;
        }

        // Query overlaps a leaf node => report collision.

        if (overlapL && lbvh.isLeaf(childL)){
          const uint32_t l = lbvh.getNumPoints(childL);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childL, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < lmin.z)
                lmin = p;
            }
          }
        }
        if (overlapR && lbvh.isLeaf(childR)){
          const uint32_t l = lbvh.getNumPoints(childR);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childR, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < lmin.z)
                lmin = p;
            }
          }
        }

        // Query overlaps an internal node => traverse.

        bool traverseL = (overlapL && !lbvh.isLeaf(childL));
        bool traverseR = (overlapR && !lbvh.isLeaf(childR));

        if (!traverseL && !traverseR)
            node = *--stackPtr; // pop
        else
        {
            node = (traverseL) ? childL : childR;
            if (traverseL && traverseR)
                *stackPtr++ = childR; // push
        }
    }
    while (node != NULL);

    // count = numPts;

    return;
}

/*A diferencia de traverseiterative2, lo que hago en este es devolver solo los índices
y no me preocupo de guardar el punto completo*/
__device__ 
void traverseIterative3(LBVH& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    uint32_t stack[64];
    uint32_t* stackPtr = stack;
    *stackPtr++ = NULL; // push

    real_t zmin = inf;

    // Traverse nodes starting from the root.

    uint32_t node = lbvh.getRoot();
    do
    {
        // Check each child node for overlap.

        uint32_t childL = lbvh.getLeftChild(node);
        uint32_t childR = lbvh.getRightChild(node);
        bool overlapL = (boxOverlap2D(queryAABB, lbvh.getAABB(childL)));
        bool overlapR = (boxOverlap2D(queryAABB, lbvh.getAABB(childR)));

        // Internal node totally overlaped

        if(isAllOverlaped(queryAABB, lbvh.getAABB(childL))){
          overlapL = false;
          numPts += lbvh.getNumPoints(childL);
          point_t min = lbvh.getMin(childL);
          if(min.z < zmin){
            zmin = min.z;
            idmin = lbvh.getMinIdx(childL);
          }
        }
        if(isAllOverlaped(queryAABB, lbvh.getAABB(childR))){
          overlapR = false;
          numPts += lbvh.getNumPoints(childR);
          point_t min = lbvh.getMin(childR);
          if(min.z < zmin){
            zmin = min.z;
            idmin = lbvh.getMinIdx(childR);
          }
        }

        // Query overlaps a leaf node => report collision.

        if (overlapL && lbvh.isLeaf(childL)){
          const uint32_t l = lbvh.getNumPoints(childL);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childL, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < zmin)
                zmin = p.z;
                idmin = lbvh.getPointIdx(childL, i);
            }
          }
        }
        if (overlapR && lbvh.isLeaf(childR)){
          const uint32_t l = lbvh.getNumPoints(childR);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childR, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < zmin)
                zmin = p.z;
                idmin = lbvh.getPointIdx(childR, i);
            }
          }
        }

        // Query overlaps an internal node => traverse.

        bool traverseL = (overlapL && !lbvh.isLeaf(childL));
        bool traverseR = (overlapR && !lbvh.isLeaf(childR));

        if (!traverseL && !traverseR)
            node = *--stackPtr; // pop
        else
        {
            node = (traverseL) ? childL : childR;
            if (traverseL && traverseR)
                *stackPtr++ = childR; // push
        }
    }
    while (node != NULL);

    // count = numPts;

    return;
}

/*A diferencia de traverseiterative3, lo que voy a hacer en este es obtener directemente
una tupla con Z e ID del punto que quiero*/
template<typename lbvh_t>
__device__ 
void traverseIterative4(lbvh_t& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    uint32_t stack[64];
    uint32_t* stackPtr = stack;
    *stackPtr++ = NULL; // push

    real_t zmin = inf;

    // Traverse nodes starting from the root.

    uint32_t node = lbvh.getRoot();
    do
    {
        // Check each child node for overlap.

        uint32_t childL = lbvh.getLeftChild(node);
        uint32_t childR = lbvh.getRightChild(node);
        bool overlapL = (boxOverlap2D(queryAABB, lbvh.getAABB(childL)));
        bool overlapR = (boxOverlap2D(queryAABB, lbvh.getAABB(childR)));

        // Internal node totally overlaped

        if(isAllOverlaped(queryAABB, lbvh.getAABB(childL))){
          overlapL = false;
          numPts += lbvh.getNumPoints(childL);
          auto min = lbvh.getMinTuple(childL);
          if(min.first < zmin){
            zmin = min.first;
            idmin = min.second;
          }
        }
        if(isAllOverlaped(queryAABB, lbvh.getAABB(childR))){
          overlapR = false;
          numPts += lbvh.getNumPoints(childR);
          auto min = lbvh.getMinTuple(childR);
          if(min.first < zmin){
            zmin = min.first;
            idmin = min.second;
          }
        }

        // Query overlaps a leaf node => report collision.

        if (overlapL && lbvh.isLeaf(childL)){
          const uint32_t l = lbvh.getNumPoints(childL);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childL, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < zmin)
                zmin = p.z;
                idmin = uint32_t(p.w);
            }
          }
        }
        if (overlapR && lbvh.isLeaf(childR)){
          const uint32_t l = lbvh.getNumPoints(childR);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childR, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < zmin)
                zmin = p.z;
                idmin = uint32_t(p.w);
            }
          }
        }

        // Query overlaps an internal node => traverse.

        bool traverseL = (overlapL && !lbvh.isLeaf(childL));
        bool traverseR = (overlapR && !lbvh.isLeaf(childR));

        if (!traverseL && !traverseR)
            node = *--stackPtr; // pop
        else
        {
            node = (traverseL) ? childL : childR;
            if (traverseL && traverseR)
                *stackPtr++ = childR; // push
        }
    }
    while (node != NULL);

    // count = numPts;

    return;
}





__global__ void query(const node_t* n, const aabb_t* bb, aabb_t initBox, double Displace, const point_t* p,
  uint32_t* count, const uint32_t* ids, uint32_t numNodes, uint32_t numBB, uint32_t minNumPoints, uint32_t nCols)
{
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx >= numBB)
      return;

  LBVH lbvh;
  lbvh.node_list       = n;
  lbvh.aabb_list       = bb;
  lbvh.objectOrder     = ids;
  lbvh.point_list      = p;
  lbvh.numNodes        = numNodes;

  aabb_t cellBox;
  cellBox.upper.x = initBox.upper.x + (idx%nCols)*Displace;
  cellBox.lower.x = initBox.lower.x + (idx%nCols)*Displace;
  cellBox.upper.y = initBox.upper.y + (int)(idx/nCols)*Displace;
  cellBox.lower.y = initBox.lower.y + (int)(idx/nCols)*Displace;

  // uint ydis = (uint)(idx%nCols);
  // uint xdis = (uint)(idx/nCols);
  // cellBox.upper.y = initBox.upper.y + ydis*Displace;
  // cellBox.lower.y = initBox.lower.y + ydis*Displace;
  // cellBox.upper.x = initBox.upper.x + xdis*Displace;
  // cellBox.lower.x = initBox.lower.x + xdis*Displace;

  uint32_t minCount = 0;
  // point_t lmin = {0.0,0.0,9999.0,0.0};
  uint32_t idmin = 0xFFFFFFFFu;

  // uint32_t min = traverseRecursive(lbvh, list[idx], lbvh.getRoot());
  // traverseIterative(lbvh, cellBox, minCount, lmin);
  // traverseIterative2(lbvh, cellBox, minCount, lmin);
  // traverseIterative3(lbvh, cellBox, minCount, idmin);
  traverseIterative4(lbvh, cellBox, minCount, idmin);

  __syncthreads();

  if( minNumPoints <= minCount){
    // count[idx] = minCount;
    count[idx] = idmin;
  }

  return;

}

__global__ void query_ordered(const node_t* n, const aabb_t* bb, aabb_t initBox, double Displace, const point_t* p,
  uint32_t* count, uint32_t Ncells, uint32_t minNumPoints, uint32_t nCols)
{
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx >= Ncells)
      return;

  // assert(idx < Ncells);

  LBVHo lbvh;
  lbvh.node_list       = n;
  lbvh.aabb_list       = bb;
  lbvh.point_list      = p;
  // lbvh.numInternalNodes        = numNodes;

  aabb_t cellBox;
  cellBox.upper.x = initBox.upper.x + (idx%nCols)*Displace;
  cellBox.lower.x = initBox.lower.x + (idx%nCols)*Displace;
  cellBox.upper.y = initBox.upper.y + (int)(idx/nCols)*Displace;
  cellBox.lower.y = initBox.lower.y + (int)(idx/nCols)*Displace;

  uint32_t minCount = 0;
  // point_t lmin = {0.0,0.0,9999.0,0.0};
  uint32_t idmin = 0xFFFFFFFFu;

  traverseIterative4(lbvh, cellBox, minCount, idmin);

  __syncthreads();

  if( minNumPoints <= minCount){
    // count[idx] = minCount;
    count[idx] = idmin;
  }

  return;

}

__global__ void query2D(const node_t* n, const aabb_t* bb, aabb_t initBox, double Displace, const point_t* p,
  uint32_t* count, const uint32_t* ids, uint32_t numNodes, uint32_t minNumPoints, uint32_t nCols, uint32_t nRows)
{
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  int jdx = threadIdx.y + blockDim.y * blockIdx.y;
  if (idx >= nCols || jdx >= nRows)
      return;

  LBVH lbvh;
  lbvh.node_list       = n;
  lbvh.aabb_list       = bb;
  lbvh.objectOrder     = ids;
  lbvh.point_list      = p;
  lbvh.numNodes        = numNodes;

  aabb_t cellBox;
  cellBox.upper.y = initBox.upper.y + jdx*Displace;
  cellBox.lower.y = initBox.lower.y + jdx*Displace;
  cellBox.upper.x = initBox.upper.x + idx*Displace;
  cellBox.lower.x = initBox.lower.x + idx*Displace;

  uint32_t minCount = 0;
  point_t lmin = {0.0,0.0,9999.0,0.0};
  // uint32_t min = traverseRecursive(lbvh, list[idx], lbvh.getRoot());
  traverseIterative(lbvh, cellBox, minCount, lmin);

  __syncthreads();

  if( minNumPoints <= minCount){
    // count[idx] = minCount;
  }

  return;

}

/* Recursive search of the minimum in an area */
point_t findValidMin(node_t* n, aabb_t* bb, point_t* p, uint32_t* idxs, 
    node_t& current, aabb_t& searchBB, uint32_t& numPts)
{
  // Lpoint tmp, min = nomin;
  point_t tmp, min = {0.0,0.0,99999.0,0.0};

  if(current.left_idx == 0xFFFFFFFFu)
  {

    const uint32_t* indices = &idxs[current.object_idx];
    uint32_t countPoints=0;
    for(int i=0; i<current.numPts; i++){
      if(insideBox(searchBB, p[indices[i]])){
        countPoints++;
        if(p[indices[i]].z - min.z < TOL)
          min = p[indices[i]];
      }
    }
    numPts += countPoints;

  } else {

    int childIdx = current.left_idx;
    if(boxOverlap2D(searchBB, bb[childIdx])){
      if(isAllOverlaped(searchBB, bb[childIdx])){
        numPts += n[childIdx].numPts;
        if (p[n[childIdx].min_idx].z - min.z < TOL) {
            min = p[n[childIdx].min_idx];
        }
      }
      else if(insideBox(searchBB,p[n[childIdx].min_idx])){
        if (p[n[childIdx].min_idx].z - min.z < TOL) {
          min = p[n[childIdx].min_idx];
        }
        countPoints(n,bb,p,idxs,n[childIdx], searchBB, numPts);
      }
      else {
        tmp = findValidMin(n,bb,p,idxs,n[childIdx], searchBB, numPts);
        if (tmp.z - min.z < TOL) {
            min = tmp;
        }
      }
    }

    childIdx = current.right_idx;
    if(boxOverlap2D(searchBB, bb[childIdx])){
      if(isAllOverlaped(searchBB, bb[childIdx])){
        numPts += n[childIdx].numPts;
        if (p[n[childIdx].min_idx].z - min.z < TOL) {
            min = p[n[childIdx].min_idx];
        }
      }
      else if(insideBox(searchBB,p[n[childIdx].min_idx])){
        if (p[n[childIdx].min_idx].z - min.z < TOL) {
          min = p[n[childIdx].min_idx];
        }
        countPoints(n,bb,p,idxs,n[childIdx], searchBB, numPts);
      }
      else {
        tmp = findValidMin(n,bb,p,idxs,n[childIdx], searchBB, numPts);
        if (tmp.z - min.z < TOL) {
            min = tmp;
        }
      }
    }

  }

  return min;
}

/*este se utiliza para recorrer el árbol con la CPU sin usar índices*/
point_t findValidMin(node_t* n, aabb_t* bb, point_t* p, 
    node_t& current, aabb_t& searchBB, uint32_t& numPts)
{
  // Lpoint tmp, min = nomin;
  point_t tmp, min = {0.0,0.0,99999.0,0.0};

  if(current.left_idx == 0xFFFFFFFFu)
  {

    point_t* points = &p[current.object_idx];
    uint32_t countPoints=0;
    for(int i=0; i<current.numPts; i++){
      if(insideBox(searchBB, points[i])){
        countPoints++;
        if(points[i].z - min.z < TOL)
          min = points[i];
      }
    }
    numPts += countPoints;

  } else {

    int childIdx = current.left_idx;
    if(boxOverlap2D(searchBB, bb[childIdx])){
      if(isAllOverlaped(searchBB, bb[childIdx])){
        numPts += n[childIdx].numPts;
        if (p[n[childIdx].min_idx].z - min.z < TOL) {
            min = p[n[childIdx].min_idx];
        }
      }
      else {
        tmp = findValidMin(n,bb,p,n[childIdx], searchBB, numPts);
        if (tmp.z - min.z < TOL) {
            min = tmp;
        }
      }
    }

    childIdx = current.right_idx;
    if(boxOverlap2D(searchBB, bb[childIdx])){
      if(isAllOverlaped(searchBB, bb[childIdx])){
        numPts += n[childIdx].numPts;
        if (p[n[childIdx].min_idx].z - min.z < TOL) {
            min = p[n[childIdx].min_idx];
        }
      }
      else {
        tmp = findValidMin(n,bb,p,n[childIdx], searchBB, numPts);
        if (tmp.z - min.z < TOL) {
            min = tmp;
        }
      }
    }

  }

  return min;
}

void checkStruct(node_t* n, node_t& current, uint32_t level, std::vector<int>& list)
{

  if(current.left_idx == 0xFFFFFFFFu)
  {

    if(list.size() <= level)
      list.push_back(1);
    else
      list[level]++;

  } else {

    if(list.size() <= level)
      list.push_back(1);
    else
      list[level]++;

    int childIdx = current.left_idx;
    checkStruct(n,n[childIdx], level+1, list);


    childIdx = current.right_idx;
    checkStruct(n,n[childIdx], level+1, list);


  }

  return;
}

template<typename pointer_t>
void checkStruct(pointer_t current, uint32_t level, std::vector<int>& list)
{

  if(current->childs == NULL)
  {

    if(list.size() <= level)
      list.push_back(1);
    else
      list[level]++;

  } else {

    if(list.size() <= level)
      list.push_back(1);
    else
      list[level]++;

    pointer_t child = &(current->childs[0]);
    checkStruct(child, level+1, list);


    child = &(current->childs[1]);
    checkStruct(child, level+1, list);


  }

  return;
}

void checkStruct(Btree current, uint32_t level, std::vector<int>& list)
{

  if(isLeaf(current))
  {

    if(list.size() <= level)
      list.push_back(1);
    else
      list[level]++;

  } else {

    if(list.size() <= level)
      list.push_back(1);
    else
      list[level]++;

    Btree child = current->childs[0];
    checkStruct(child, level+1, list);


    child = current->childs[1];
    checkStruct(child, level+1, list);


  }

  return;
}

void checkArea(tree2 current, std::vector<std::pair<Vector2D,float>>& areas)
{

  if(current->childs == NULL)
  {

    areas.push_back(std::make_pair(current->center, current->radius));

  } else {

    tree2 child = &(current->childs[0]);
    checkArea(child, areas);


    child = &(current->childs[1]);
    checkArea(child, areas);


  }

  return;
}


/* Stage1 with parallel_for that remembers previous minimum */
void stage1tbbRem(node_t* n, aabb_t* bb, point_t* p, uint32_t* idxs,
  uint32_t Wsize, double Overlap, uint32_t nCols, uint32_t nRows,
  uint32_t minNumPoints, point_t* minVector, node_t& root, aabb_t& BBox)
{

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = BBox.lower.x - Wsize + Displace;
  initBox.lower.y = BBox.lower.y - Wsize + Displace;
  initBox.upper.x = BBox.lower.x + Displace;
  initBox.upper.y = BBox.lower.y + Displace;

  // tbb::blocked_range<int>::size_type gs = 4;

  // tbb::affinity_partitioner aff_p;

  tbb::parallel_for( tbb::blocked_range2d<uint,uint>{0,nRows,0,nCols},
                      [&](tbb::blocked_range2d<uint,uint> r ) {

    point_t newmin = {0,0.0,0.0,0.0};
    // Vector2D cellCenter;
    aabb_t cellBox;
    // int cellPoints = 0;
    int je = r.rows().end();
    int ie = r.cols().end();

    for(int jj = r.rows().begin(); jj < je; ++jj) {

      uint32_t cellPoints = 0;

      cellBox.upper.y = initBox.upper.y + jj*Displace;
      cellBox.lower.y = initBox.lower.y + jj*Displace;

      for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

        cellBox.upper.x = initBox.upper.x + ii*Displace;
        cellBox.lower.x = initBox.lower.x + ii*Displace;

        // makeBox(cellCenter, Wsize*0.5, boxMin, boxMax);

        if(insideBox(cellBox, newmin)){
        // if(cellPoints > 0 && insideBox(&newmin,boxMin,boxMax)){

          // Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};
          aabb_t overlapBox;
          overlapBox.upper.x = cellBox.upper.x;
          overlapBox.upper.y = cellBox.upper.y;
          overlapBox.lower.x = cellBox.upper.x - Displace;
          overlapBox.lower.y = cellBox.lower.y;

          int old_cellPoints = cellPoints;
          cellPoints = 0;

          point_t tmp = findValidMin(n, bb, p, idxs, root, overlapBox, cellPoints);
          // point_t tmp = findValidMin(n, bb, p, root, overlapBox, cellPoints);

          // We're assuming the points were equidistant throughout the cell, which isn't always true.
          cellPoints += (int)(old_cellPoints * Overlap);

          if(tmp.z < newmin.z){
            newmin = tmp;
          }

        } else {

          cellPoints = 0;
          newmin = findValidMin(n, bb, p, idxs, root, cellBox, cellPoints);
          // newmin = findValidMin(n, bb, p, root, cellBox, cellPoints);

        }

        if( cellPoints >= minNumPoints ){

            // v.push_back(newmin.id);
            minVector[jj*nCols + ii] = newmin;

        }

      }
    }

  });

  return;
}

/* Stage1 with parallel_for that remembers previous minimum */
void stage1tbb(node_t* n, aabb_t* bb, point_t* p, uint32_t* idxs, uint32_t* count,
  uint32_t Wsize, double Overlap, uint32_t nCols, uint32_t nRows,
  uint32_t minNumPoints, aabb_t& BBox, uint32_t numThreads, uint32_t numNodes)
{

  size_t Ncells = nRows*nCols;

  cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = BBox.lower.x - Wsize + Displace;
  initBox.lower.y = BBox.lower.y - Wsize + Displace;
  initBox.upper.x = BBox.lower.x + Displace;
  initBox.upper.y = BBox.lower.y + Displace;

  int numBlocks = int((Ncells-1)/numThreads) + 1;

  query<<<numBlocks, numThreads>>>(n, bb, initBox, Displace, p, count, idxs, numNodes, Ncells, minNumPoints, nCols);

  // dim3 tblocks(16,16,1);
  // dim3 grid((nCols/tblocks.x)+1, (nRows/tblocks.y)+1, 1);
  // query2D<<<grid, tblocks>>>(n, bb, initBox, Displace, p, count, idxs, numNodes, minNumPoints, nCols, nRows);

  cudaDeviceSynchronize();

  // cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);

  // for(int i=0; i<Ncells; i++){
  //   if(count[i] != 0)
  //     minVector[i].x = 1.0;
  // }

  // cudaFree(list);
  // cudaFree(count);

  return;
}

void stage1tbb(node_t* n, aabb_t* bb, point_t* p, uint32_t* count,
  uint32_t Wsize, double Overlap, uint32_t nCols, uint32_t nRows,
  uint32_t minNumPoints, aabb_t& BBox, uint32_t numThreads)
{

  size_t Ncells = nRows*nCols;

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = BBox.lower.x - Wsize + Displace;
  initBox.lower.y = BBox.lower.y - Wsize + Displace;
  initBox.upper.x = BBox.lower.x + Displace;
  initBox.upper.y = BBox.lower.y + Displace;

  int numBlocks = int((Ncells-1)/numThreads) + 1;

  query_ordered<<<numBlocks, numThreads>>>(n, bb, initBox, Displace, p, count, Ncells, minNumPoints, nCols);

  cudaDeviceSynchronize();


  return;
}


/* Recursively copy the tree using pointers */
void copyTree(const node_t* n, const aabb_t* bb, point_t* p, const uint32_t* idxs, 
  tree bt, const node_t& current)
{
  if(current.left_idx == 0xFFFFFFFFu) {

    const uint vsize = current.numPts;

    uint index = num_pts_copiados;

    bt->points = &(array_all_points[num_pts_copiados]);

    bt->numPts = vsize;

    bt->min = &(p[current.min_idx]);

    num_pts_copiados += vsize;

    const uint32_t* indices = &idxs[current.object_idx];
    for(int i=0; i<current.numPts; i++){
      array_all_points[index] = p[indices[i]];
      index++;
    }



  } else {

    // ahora es cuando se que quadrants debe apuntar a algo
    bt->childs = &(array_nodes[copied_nodes]);

    // copio el número de puntos
    bt->numPts = current.numPts;

    /* copio el minimo; TODO: esto está mal porque debería ser el mínimo en el nuevo vector, 
    pero el problema es que se el mínimo, pero no su posición en el nuevo vector. Debería
    dejar la comprobación del mínimo para esta etapa o hacerlo de forma diferente*/
    bt->min = &(p[current.min_idx]);

    tree dst = &(bt->childs[0]);

    new (dst) tree_t( bt, bb[current.left_idx]);

    dst = &(bt->childs[1]);

    new (dst) tree_t( bt, bb[current.right_idx]);

    copied_nodes += 2;


    // dst = &(bt->childs[0]);

    copyTree(n, bb, p, idxs, &(bt->childs[0]), n[current.left_idx]);

    // dst = &(bt->childs[1]);

    copyTree(n, bb, p, idxs, &(bt->childs[1]), n[current.right_idx]);

  }

}

tree launchCopy(const node_t* n, const aabb_t* bb, point_t* p, const uint32_t* idxs)
{
  tree root = &(array_nodes[0]);

  new (root) tree_t(  NULL, bb[0] );

  copied_nodes++; // hasta aquí ya tengo copiado root

  copyTree(n, bb, p, idxs, root, n[0]);
  
  return root;
  
}

/* Recursively copy the tree using pointers */
void copyTree2(const node_t* n, const aabb_t* bb, point_t* p, const uint32_t* idxs, 
  tree bt, const node_t& current)
{
  if(current.left_idx == 0xFFFFFFFFu) {

    const uint vsize = current.numPts;

    uint index = num_pts_copiados;

    bt->points = &(array_all_points[num_pts_copiados]);

    bt->numPts = vsize;

    // bt->min = &(p[current.min_idx]);

    num_pts_copiados += vsize;

    point_t* min = NULL;
    const uint32_t* indices = &idxs[current.object_idx];
    for(int i=0; i<current.numPts; i++){
      array_all_points[index] = p[indices[i]];

      if(min == NULL) {
        // min = p;
        min = &(array_all_points[index]);
      } else if(p[indices[i]].z < min->z) {
        min = &(array_all_points[index]);
      }
    
      index++;
    }

    bt->min = min;


  } else {

    // ahora es cuando se que quadrants debe apuntar a algo
    bt->childs = &(array_nodes[copied_nodes]);

    // copio el número de puntos
    bt->numPts = current.numPts;

    /* copio el minimo; TODO: esto está mal porque debería ser el mínimo en el nuevo vector, 
    pero el problema es que se el mínimo, pero no su posición en el nuevo vector. Debería
    dejar la comprobación del mínimo para esta etapa o hacerlo de forma diferente*/
    // bt->min = &(p[current.min_idx]);

    tree dst = &(bt->childs[0]);

    new (dst) tree_t( bt, bb[current.left_idx]);

    dst = &(bt->childs[1]);

    new (dst) tree_t( bt, bb[current.right_idx]);

    copied_nodes += 2;


    point_t* min = NULL;
    dst = &(bt->childs[0]);

    copyTree2(n, bb, p, idxs, dst, n[current.left_idx]);

    if(dst->min != NULL) {
      if(min == NULL) {
        min = dst->min;
      } else if(dst->min->z < min->z) {
        min = dst->min;
      }
    }

    dst = &(bt->childs[1]);

    copyTree2(n, bb, p, idxs, dst, n[current.right_idx]);

    if(dst->min != NULL) {
      if(min == NULL) {
        min = dst->min;
      } else if(dst->min->z < min->z) {
        min = dst->min;
      }
    }

    bt->min = min;

  }

}

tree launchCopy2(const node_t* n, const aabb_t* bb, point_t* p, const uint32_t* idxs)
{
  tree root = &(array_nodes[0]);

  new (root) tree_t(  NULL, bb[0] );

  copied_nodes++; // hasta aquí ya tengo copiado root

  copyTree2(n, bb, p, idxs, root, n[0]);
  
  return root;
  
}



node_t* array_nodes_gpu = NULL;
aabb_t* array_aabbs_gpu = NULL;
uint32_t gpu_copied_nodes = 0;

bool checkNodes(const node_t* n, const aabb_t* bb, point_t* p, uint32_t current_idx, uint32_t oldCurrent_idx)
{
  node_t* current = &(array_nodes_gpu[current_idx]);

  const node_t* oldCurrent = &(n[oldCurrent_idx]);

  if( current->numPts != oldCurrent->numPts || current->object_idx != oldCurrent->object_idx ){
        printf("fallo en el nodo\n");
        return 0;
  }


  if( fabs(array_aabbs_gpu[current_idx].upper.x - bb[oldCurrent_idx].upper.x) > TOL ||
      fabs(array_aabbs_gpu[current_idx].upper.y - bb[oldCurrent_idx].upper.y) > TOL ||
      fabs(array_aabbs_gpu[current_idx].lower.x - bb[oldCurrent_idx].lower.x) > TOL ||
      fabs(array_aabbs_gpu[current_idx].lower.y - bb[oldCurrent_idx].lower.y) > TOL ){
        printf("fallo BB\n");
        return 0;
      }


  if( /*array_all_points[current->min_idx].x != p[oldCurrent->min_idx].x ||
      array_all_points[current->min_idx].y != p[oldCurrent->min_idx].y ||*/
      array_all_points[current->min_idx].z != p[oldCurrent->min_idx].z ){
        printf("fallo en el minimo\n");
        printf("current_idx: %d, oldCurrent_idx: %d\n", current_idx, oldCurrent_idx);
        printf("current_min_idx: %d, oldCurrent_min_idx: %d\n", current->min_idx, oldCurrent->min_idx);
        printf("current min: %lf, %lf, %lf\n", array_all_points[current->min_idx].x, array_all_points[current->min_idx].y, array_all_points[current->min_idx].z);
        printf("oldCurrent min: %lf, %lf, %lf\n", p[oldCurrent->min_idx].x, p[oldCurrent->min_idx].y, p[oldCurrent->min_idx].z);
        return 0;
      }

  if(current->min_idx >= num_objects){
    printf("fallo en el indice\n");
    return 0;
  }
  if(current_idx!=0)
    if(current->parent_idx >= r_num_nodes ||
      current->left_idx >= r_num_nodes   ||
      current->right_idx >= r_num_nodes ){
      printf("fallo indices nodos\n");
      printf("num_nodes = %u\n", r_num_nodes);
      printf("parent_idx = %u\n", current->parent_idx);
      printf("left_idx = %u\n", current->left_idx);
      printf("right_idx = %u\n", current->right_idx);
      return 0;
    }

  return 1;
}

void copyGPU(const node_t* n, const aabb_t* bb, point_t* p, const uint32_t* idxs, uint32_t current_idx, uint32_t oldCurrent_idx)
{
  node_t* current = &(array_nodes_gpu[current_idx]);

  const node_t* oldCurrent = &(n[oldCurrent_idx]);

  if(oldCurrent->left_idx == 0xFFFFFFFFu) { //isLeaf??

    uint32_t length = oldCurrent->numPts;

    uint32_t index = num_pts_copiados;

    // establezco la posición al primer punto
    current->object_idx = num_pts_copiados;

    // copio el numero de puntos
    current->numPts = length;

    // copio la BB
    array_aabbs_gpu[current_idx] = bb[oldCurrent_idx];

    // reservo el número de puntos que voy a copiar
    num_pts_copiados += length;

    uint32_t idmin = 0xFFFFFFFFu;
    const uint32_t* indices = &idxs[oldCurrent->object_idx];

    // copio los puntos y busco el mínimo
    for(int i=0; i<length; i++){

      array_all_points[index] = p[indices[i]];

      if(idmin == 0xFFFFFFFFu) {
        idmin = index;
      } else if(p[indices[i]].z < array_all_points[idmin].z) {
        idmin = index;
      }
    
      index++;
    }

    // ya conozco el minimo
    assert(idmin != 0xFFFFFFFFu);
    current->min_idx = idmin;


  } else {

    uint32_t idmin = 0xFFFFFFFFu;
    // ahora es cuando se que childs debe apuntar a algo
    current->left_idx = gpu_copied_nodes;
    current->right_idx = gpu_copied_nodes+1;

    gpu_copied_nodes += 2;

    // copio el número de puntos
    current->numPts = oldCurrent->numPts;

    // copio la BB
    array_aabbs_gpu[current_idx] = bb[oldCurrent_idx];

    // copio el primer nodo hijo
    node_t* dst = &(array_nodes_gpu[current->left_idx]);

    *(dst) = default_node;

    dst->parent_idx = current_idx;

    copyGPU(n, bb, p, idxs, current->left_idx, oldCurrent->left_idx);

    if(dst->min_idx != 0xFFFFFFFFu)
    {
      if(idmin == 0xFFFFFFFFu){
        idmin = dst->min_idx;
      } else if(array_all_points[dst->min_idx].z < array_all_points[idmin].z){
        idmin = dst->min_idx;
      }
    }

    // copio el segundo
    dst = &(array_nodes_gpu[current->right_idx]);

    *(dst) = default_node;

    dst->parent_idx = current_idx;

    copyGPU(n, bb, p, idxs, current->right_idx, oldCurrent->right_idx);

    if(dst->min_idx != 0xFFFFFFFFu)
    {
      if(idmin == 0xFFFFFFFFu){
        idmin = dst->min_idx;
      } else if(array_all_points[dst->min_idx].z < array_all_points[idmin].z){
        idmin = dst->min_idx;
      }
    }

    // ya conozco el mínimo
    assert(idmin != 0xFFFFFFFFu);
    current->min_idx = idmin;

    // assert(checkNodes(n, bb, p, current_idx, oldCurrent_idx));
  
  }
  
}

void launchCopyGPU(const node_t* n, const aabb_t* bb, point_t* p, const uint32_t* idxs)
{

  array_nodes_gpu[0] = default_node;

  gpu_copied_nodes++; // hasta aquí ya tengo copiado root

  copyGPU(n, bb, p, idxs, 0, 0);

  return;
}

float getRadius(const aabb_t& bb){
  
  float radioX = (bb.upper.x - bb.lower.x)*0.5;
  float radioY = (bb.upper.y - bb.lower.y)*0.5;

  return (radioX < radioY)? radioY : radioX;
}

Vector2D getCenter(const aabb_t& bb)
{
    Vector2D center;

    center.x = bb.lower.x + (bb.upper.x - bb.lower.x)*0.5;
    center.y = bb.lower.y + (bb.upper.y - bb.lower.y)*0.5;
    // center.z = min.z + radius.z;

    return center;
}

/* Recursively copy the tree using pointers */
template<typename pointer_t>
void copyTree3(const node_t* n, const aabb_t* bb, point_t* p, const uint32_t* idxs, 
  pointer_t bt, const node_t& current)
{
  if(current.left_idx == 0xFFFFFFFFu) {

    const uint vsize = current.numPts;

    uint index = num_pts_copiados;

    bt->points = &(array_all_points[num_pts_copiados]);

    bt->numPts = vsize;

    // bt->min = &(p[current.min_idx]);

    num_pts_copiados += vsize;

    point_t* min = NULL;
    const uint32_t* indices = &idxs[current.object_idx];
    for(int i=0; i<current.numPts; i++){
      array_all_points[index] = p[indices[i]];

      if(min == NULL) {
        // min = p;
        min = &(array_all_points[index]);
      } else if(p[indices[i]].z < min->z) {
        min = &(array_all_points[index]);
      }
    
      index++;
    }

    bt->min = min;


  } else {

    // ahora es cuando se que quadrants debe apuntar a algo
    bt->childs = &(array_nodes2[copied_nodes]);

    // copio el número de puntos
    bt->numPts = current.numPts;

    /* copio el minimo; TODO: esto está mal porque debería ser el mínimo en el nuevo vector, 
    pero el problema es que se el mínimo, pero no su posición en el nuevo vector. Debería
    dejar la comprobación del mínimo para esta etapa o hacerlo de forma diferente*/
    // bt->min = &(p[current.min_idx]);

    pointer_t dst = &(bt->childs[0]);

    new (dst) tree2_t( bt, getCenter(bb[current.left_idx]), getRadius(bb[current.left_idx]) );

    dst = &(bt->childs[1]);

    new (dst) tree2_t( bt, getCenter(bb[current.right_idx]), getRadius(bb[current.right_idx]) );

    copied_nodes += 2;


    point_t* min = NULL;
    dst = &(bt->childs[0]);

    copyTree3(n, bb, p, idxs, dst, n[current.left_idx]);

    if(dst->min != NULL) {
      if(min == NULL) {
        min = dst->min;
      } else if(dst->min->z < min->z) {
        min = dst->min;
      }
    }

    dst = &(bt->childs[1]);

    copyTree3(n, bb, p, idxs, dst, n[current.right_idx]);

    if(dst->min != NULL) {
      if(min == NULL) {
        min = dst->min;
      } else if(dst->min->z < min->z) {
        min = dst->min;
      }
    }

    bt->min = min;

  }

}


tree2 launchCopy3(const node_t* n, const aabb_t* bb, point_t* p, const uint32_t* idxs)
{
  tree2 root = &(array_nodes2[0]);

  new (root) tree2_t(  NULL, getCenter(bb[0]), getRadius(bb[0]) );

  copied_nodes++; // hasta aquí ya tengo copiado root

  copyTree3(n, bb, p, idxs, root, n[0]);
  
  return root;
  
}

/* Recursively copy the tree using pointers */
template<typename cpunode_t>
void copyTree4(cpunode_t bt, tree2 parent)
{
  if(isLeaf(bt)) {

    // si este nodo es hoja, solo tengo que copiar los puntos ya que lo he guardado anteriormente
    // quadrants lo dejo a NULL porque es hoja; evito ocupar más espacio

    size_t vsize = bt->points.size();

    if(vsize > 0) {

      uint64_t index = num_pts_copiados;

      parent->points = &(array_all_points[num_pts_copiados]);
      
      parent->numPts = vsize;

      num_pts_copiados += vsize;

      point_t* min = NULL;
      for(point_t* p : bt->points) {
        array_all_points[index] = *p;
      
        if(min == NULL) {
          // min = p;
          min = &(array_all_points[index]);
        } else if(p->z < min->z) {
          min = &(array_all_points[index]);
        }

        index++;
      }

      parent->min = min;
    }

  } else {

    // ahora es cuando se que quadrants debe apuntar a algo
    parent->childs = &(array_nodes2[copied_nodes]);

    // copio los 4 cuadrantes contiguamente
    for(int i=0; i<2; i++) {

      cpunode_t src = bt->childs[i];

      tree2 dst = &(parent->childs[i]);

      new (dst) tree2_t(  parent, 
                            src->center, 
                            src->radius );
      
      // copied_nodes++;
    }

    copied_nodes += 2;

    // continúo recursivamente para cada uno de los cuadrantes
    point_t* min = NULL;
    for(int i=0; i<2; i++) {

      tree2 dst = &(parent->childs[i]);

      copyTree4(bt->childs[i], dst);

      parent->numPts += dst->numPts;

      if(dst->min != NULL) {
        if(min == NULL) {
          min = dst->min;
        } else if(dst->min->z < min->z) {
          min = dst->min;
        }
      }

    }

    parent->min = min;

  }

}

template<typename cpunode_t>
tree2 launchCopy4(cpunode_t bt)
{
  tree2 root = &(array_nodes2[0]);

  new (root) tree2_t(  NULL, 
                        bt->center, 
                        bt->radius );

  copied_nodes++; // hasta aquí ya tengo copiado root

  copyTree4(bt, root);
  
  return root;
  
}

/* Recursive search of the minimum in an area */
void countPoints(tree current, aabb_t& searchBB, uint32_t& numPts)
{

  if(current->childs == NULL)
  {
    uint32_t size = current->numPts;
    point_t* p = current->points;
    uint32_t countPoints=0;
    for(int i=0; i<size; i++){
      if(insideBox(searchBB, p[i])){
        countPoints++;
      }
    }
    numPts += countPoints;

  } else {

    for(int i=0; i<2; i++)
    {
      tree child = &(current->childs[i]);
      if(boxOverlap2D(searchBB, child->bbox)){
        if(isAllOverlaped(searchBB, child->bbox)){
          numPts += child->numPts;
        }
        else {
          countPoints(child, searchBB, numPts);
        }
      }
    }

  }

  return;
}

/* Recursive search of the minimum in an area */
point_t findValidMin(tree current, aabb_t& searchBB, uint32_t& numPts)
{
  // Lpoint tmp, min = nomin;
  point_t tmp, min = {0.0,0.0,99999.0,0.0};

  if(current->childs == NULL)
  {

    uint32_t size = current->numPts;
    point_t* p = current->points;
    uint32_t countPoints=0;
    for(int i=0; i<size; i++){
      if(insideBox(searchBB, p[i])){
        countPoints++;
        if(p[i].z - min.z < TOL)
          min = p[i];
      }
    }
    numPts += countPoints;

  } else {

    for(int i=0; i<2; i++)
    {
      tree child = &(current->childs[i]);
      if(boxOverlap2D(searchBB, child->bbox)){
        if(isAllOverlaped(searchBB, child->bbox)){
          numPts += child->numPts;
          if (child->min->z - min.z < TOL) {
              min = *(child->min);
          }
        }
        else if(insideBox(searchBB,*(child->min))){
          if (child->min->z - min.z < TOL) {
            min = *(child->min);
          }
          countPoints(child, searchBB, numPts);
        }
        else {
          tmp = findValidMin(child, searchBB, numPts);
          if (tmp.z - min.z < TOL) {
              min = tmp;
          }
        }
      }
    }

  }

  return min;
}

template<typename pointer_t>
inline void findPosition(pointer_t parent, pointer_t& current, int& idx) 
{

  if(parent->childs == current){
    // std::advance(iterador, 1);
    idx=1;
    current = &(parent->childs[1]);
  }
  else{
    // std::advance(iterador, 2);
    idx=2;
    current = parent->childs;
  }

}


/* Performs the search when the quadrants are contiguos in memory */
point_t gpuSearchNeighborsMin(tree root, aabb_t& searchBB, uint& numNeighs)
{
    int idx = 0;
    uint numInside = 0;

    point_t min = {0.0,0.0,99999.0,0.0};

    tree parent = root;
    tree current = root->childs;

    while(current != root) {


      if(idx > 1) {

        current = current->parent;

        if(current != root){

          parent = current->parent;

          findPosition(parent, current, idx);

        }
        
      } else {

        if( current->childs == NULL ) { // isLeaf??

          if( current->numPts > 0 ) { //isEmpty??

            /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
            if(isAllOverlaped(searchBB, current->bbox)) {

              if(current->min->z < min.z)
                min = *(current->min);
              numInside += current->numPts;

            /*Si está parcialmente solapada, tengo dos opciones*/
            } else if(boxOverlap2D(searchBB, current->bbox)) {

              int N = current->numPts;
              point_t* pointer = current->points;

              /*Si el mínimo está dentro solo comparo una vez y cuento */
              if( insideBox(searchBB, *(current->min)) ) {

                if(current->min->z < min.z)
                  min = *(current->min);

                for(int i=0; i<N; ++i) {
                  if(insideBox(searchBB, pointer[i]))
                  {
                    numInside++;
                  }
                }
              }
              /*Si el mínimo está fuera, tengo que comprobarlo todo*/
              else {
                for(int i=0; i<N; ++i) {
                  if(insideBox(searchBB, pointer[i]))
                  {
                    if (pointer[i].z < min.z) {
                        min = pointer[i];
                    }
                    numInside++;
                  }
                }
              }

            }
            // else {
            //   std::cout << "Ni lo uno ni lo otro\n" << std::flush;
            // }

          }

          idx++;
          if(idx < 2) {
            current = &(parent->childs[idx]);

          }

        } else {
          
          if(!boxOverlap2D(searchBB, current->bbox) /* || current->numPts == 0*/ ) { //No solapada o vacia

            idx++;
            if(idx < 2) {
              current = &(parent->childs[idx]);

            }
          }
          /*si la caja está completamente solapada, capturo el mínmo
           y el numero de puntos, y continuo*/
          else if( isAllOverlaped(searchBB, current->bbox) ) {

            if(current->min->z < min.z)
              min = *(current->min);
            numInside += current->numPts;

            idx++;
            if(idx < 2) {
              current = &(parent->childs[idx]);

            }
          }
          /*pero si solo está parcialmente solapada, profundizo*/
          else {

            idx = 0;
            parent = current;
            current = current->childs;

          }
        }
      }

    }

    numNeighs = numInside;

    return min;
}


/* Recursive search of the minimum in an area */
point_t findValidMinSimple( tree current, aabb_t searchBB, uint32_t& numPts)
{
  // Lpoint tmp, min = nomin;
  point_t tmp, min = {0.0,0.0,99999.0,0.0};

  if(current->childs == NULL)
  {

  } else {

    for(int i=0; i<2; i++){
      
      tree child = &(current->childs[i]);
      if(boxOverlap2D(searchBB, child->bbox)){

        tmp = findValidMinSimple(child, searchBB, numPts);

      }
    }

  }

  return min;
}


/* Stage1 with parallel_for that remembers previous minimum */
void stage1tbbRem(uint32_t Wsize, double Overlap, uint32_t nCols, uint32_t nRows,
  uint32_t minNumPoints, point_t* minVector, tree root, aabb_t& BBox)
{

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = BBox.lower.x - Wsize + Displace;
  initBox.lower.y = BBox.lower.y - Wsize + Displace;
  initBox.upper.x = BBox.lower.x + Displace;
  initBox.upper.y = BBox.lower.y + Displace;

  // tbb::blocked_range<int>::size_type gs = 4;

  // tbb::affinity_partitioner aff_p;

  tbb::parallel_for( tbb::blocked_range2d<uint,uint>{0,nRows,0,nCols},
                      [&](tbb::blocked_range2d<uint,uint> r ) {

    point_t newmin = {0,0.0,99999.0,0.0};
    // Vector2D cellCenter;
    aabb_t cellBox;
    // int cellPoints = 0;
    int je = r.rows().end();
    int ie = r.cols().end();

    for(int jj = r.rows().begin(); jj < je; ++jj) {

      uint32_t cellPoints = 0;

      cellBox.upper.y = initBox.upper.y + jj*Displace;
      cellBox.lower.y = initBox.lower.y + jj*Displace;

      for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

        cellBox.upper.x = initBox.upper.x + ii*Displace;
        cellBox.lower.x = initBox.lower.x + ii*Displace;

        // makeBox(cellCenter, Wsize*0.5, boxMin, boxMax);

        if(insideBox(cellBox, newmin)){
        // if(cellPoints > 0 && insideBox(&newmin,boxMin,boxMax)){

          // Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};
          aabb_t overlapBox;
          overlapBox.upper.x = cellBox.upper.x;
          overlapBox.upper.y = cellBox.upper.y;
          overlapBox.lower.x = cellBox.upper.x - Displace;
          overlapBox.lower.y = cellBox.lower.y;

          int old_cellPoints = cellPoints;
          cellPoints = 0;

          point_t tmp = findValidMin(root, overlapBox, cellPoints);
          // point_t tmp = findValidMinSimple(root, overlapBox, cellPoints);
          // point_t tmp = gpuSearchNeighborsMin(root, overlapBox, cellPoints);
          // point_t tmp = {0,0.0,0.0,0.0};

          // We're assuming the points were equidistant throughout the cell, which isn't always true.
          cellPoints += (int)(old_cellPoints * Overlap);

          if(tmp.z - newmin.z < TOL){
            newmin = tmp;
          }

        } else {

          cellPoints = 0;
          newmin = findValidMin(root, cellBox, cellPoints);
          // newmin = findValidMinSimple(root, cellBox, cellPoints);
          // newmin = gpuSearchNeighborsMin(root, cellBox, cellPoints);
          // newmin = {0,0.0,0.0,0.0};

        }

        if( cellPoints >= minNumPoints ){

            // v.push_back(newmin.id);
            minVector[jj*nCols + ii] = newmin;

        }

      }
    }

  });

  return;
}

/* Makes a box of equal sides with two points */
void makeBox(Vector2D& center, float radius, Vector2D& min, Vector2D& max)
{
    min.x = center.x - radius;
    min.y = center.y - radius;
    // min.z = center.z - radius;

    max.x = center.x + radius;
    max.y = center.y + radius;
    // max.z = point.z + radius;
}

/* Makes a box of different sides with two points */
void makeBox(Vector2D& center, double radiusX, double radiusY, Vector2D& min, Vector2D& max)
{

    min.x = center.x - radiusX;
    min.y = center.y - radiusY;

    max.x = center.x + radiusX;
    max.y = center.y + radiusY;

}

/* Checks if a point is contained in the search area */
int insideBox2D(point_t* point, Vector2D& min, Vector2D& max)
{
    if(point->x > min.x && point->y > min.y)
    {
        if(point->x < max.x && point->y < max.y)
        {
            return 1;
        }
    }

    return 0;
}

/* Checks if the node is partially overlapped by the search area */
template<typename pointer_t>
int boxOverlap2D(Vector2D& boxMin, Vector2D& boxMax, pointer_t qt)
{
    if(qt->center.x + qt->radius < boxMin.x ||
       qt->center.y + qt->radius < boxMin.y)
        return 0;

    if(qt->center.x - qt->radius > boxMax.x ||
       qt->center.y - qt->radius > boxMax.y)
        return 0;

    return 1;
}


/* Checks if the search area is completely overlapped by the node */
template<typename pointer_t>
int boxInside2D(Vector2D& boxMin, Vector2D& boxMax, pointer_t qt)
{
    if(qt->center.x + qt->radius > boxMax.x ||
       qt->center.y + qt->radius > boxMax.y)
        return 0;

    if(qt->center.x - qt->radius < boxMin.x ||
       qt->center.y - qt->radius < boxMin.y)
        return 0;

    return 1;
}


/* Performs the search when the quadrants are contiguos in memory */
point_t gpuSearchNeighborsMin(Vector2D& center, tree2 btree, float radiusX, float radiusY, int& numNeighs)
{
    int idx = 0;
    int numInside = 0;
    Vector2D boxMin, boxMax;

    // *numNeighs = 0;
    makeBox(center, radiusX, radiusY, boxMin, boxMax);

    point_t min = {0.0,0.0,99999.0,0.0};

    tree2 parent = btree;
    tree2 current = btree->childs;

    while(current != btree) {


      if(idx > 1) {

        current = current->parent;

        if(current != btree){

          parent = current->parent;

          findPosition(parent, current, idx);

        }
        
      } else {

        if( current->childs == NULL ) { // isLeaf??

          if( current->numPts > 0 ) { //isEmpty??

            /*Si está completamente solapada, solo tengo que comprobar 1 punto*/
            if(boxInside2D(boxMin, boxMax, current)) {

              if(current->min->z < min.z)
                min = *(current->min);
              numInside += current->numPts;

            /*Si está parcialmente solapada, tengo dos opciones*/
            } else if(boxOverlap2D(boxMin, boxMax, current)) {

              int N = current->numPts;
              point_t* pointer = current->points;

              /*Si el mínimo está dentro solo comparo una vez y cuento */
              if( insideBox2D(current->min, boxMin, boxMax) ) {

                if(current->min->z < min.z)
                  min = *(current->min);

                for(int i=0; i<N; ++i) {
                  if(insideBox2D(&pointer[i], boxMin, boxMax))
                  {
                    numInside++;
                  }
                }
              }
              /*Si el mínimo está fuera, tengo que comprobarlo todo*/
              else {
                for(int i=0; i<N; ++i) {
                  if(insideBox2D(&pointer[i], boxMin, boxMax))
                  {
                    if (pointer[i].z < min.z) {
                        min = pointer[i];
                    }
                    numInside++;
                  }
                }
              }

            }
            // else {
            //   std::cout << "Ni lo uno ni lo otro\n" << std::flush;
            // }

          }

          idx++;
          if(idx < 2) {
            current = &(parent->childs[idx]);

          }

        } else {
          
          if(!boxOverlap2D(boxMin, boxMax, current) || current->numPts == 0) { //No solapada o vacia

            idx++;
            if(idx < 2) {
              current = &(parent->childs[idx]);

            }
          }
          /*si la caja está completamente solapada, capturo el mínmo
           y el numero de puntos, y continuo*/
          else if( boxInside2D(boxMin, boxMax, current) ) {

            if(current->min->z < min.z)
              min = *(current->min);
            numInside += current->numPts;

            idx++;
            if(idx < 2) {
              current = &(parent->childs[idx]);

            }
          }
          /*pero si solo está parcialmente solapada, profundizo*/
          else {

            idx = 0;
            parent = current;
            current = current->childs;

          }
        }
      }

    }

    numNeighs = numInside;

    return min;
}



/* Stage1 with parallel_for that remembers previous minimum */
void stage1tbbRem2(unsigned short Wsize, double Overlap, unsigned short nCols, unsigned short nRows,
  unsigned short minNumPoints, point_t* minIDs, tree2 btreeIn, double2 min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // tbb::blocked_range<int>::size_type gs = 4;

    // tbb::affinity_partitioner aff_p;

    tbb::parallel_for( tbb::blocked_range2d<int,int>{0,nRows,0,nCols},
                       [&](tbb::blocked_range2d<int,int> r ) {

        point_t newmin = {0.0,0.0,9999.0,0.0};
        Vector2D cellCenter;
        Vector2D boxMax, boxMin;
        int cellPoints = 0;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

              cellCenter.x = initX + ii*Displace;

              makeBox(cellCenter, Wsize*0.5, boxMin, boxMax);

              if(insideBox2D(&newmin,boxMin,boxMax)){
              // if(cellPoints > 0 && insideBox2D(&newmin,boxMin,boxMax)){

                Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

                int old_cellPoints = cellPoints;

                point_t tmp = gpuSearchNeighborsMin(oCell, btreeIn, Displace*0.5, Wsize*0.5, cellPoints);

                // We're assuming the points were equidistant throughout the cell, which isn't always true.
                cellPoints += (int)(old_cellPoints * Overlap);

                if(tmp.z < newmin.z){
                  newmin = tmp;
                }


              } else {

                // newmin = searchNeighborsMin(&cellCenter, btreeIn, Wsize/2, &cellPoints);
                newmin = gpuSearchNeighborsMin(cellCenter, btreeIn, Wsize*0.5, Wsize*0.5, cellPoints);

              }

              if( cellPoints >= minNumPoints ){

                  // v.push_back(newmin.id);
                  minIDs[jj*nCols + ii] = newmin;

              }

            }
        }

    });

    return;
}

template<typename pointer_t>
int isLeaf(pointer_t bt)
{
    return bt->childs[0] == NULL;
}

template<typename pointer_t>
int childIdx(point_t *point, pointer_t btree)
{
    int child = 0;

    if(point->x >= btree->center.x) child |= 1;
    // if(point->y >= btree->center.y) child |= 1;
    // if(point->z >= btree->center.z) child |= 1;

    return child;
}

uint32_t cpu_tree_nodes = 1;

void insertPoint2(point_t *point, Btree btree, int maxSize);

// /* Re-inserts the points when its exceed the set max. size  */
void fillChilds(Btree btree, int maxSize)
{

    for(point_t* p : btree->points)
    {
      int idx = childIdx(p, btree);
      insertPoint2(p, btree->childs[idx], maxSize);
    }

    btree->points.clear();
}

// /* Initializes the four quadrants */
void createChilds(Btree bt)
{
    Vector2D newCenter;
    float newRadius = bt->radius * 0.5;

    for( int i = 0; i < 2; i++)
    {
      newCenter = bt->center;
      newCenter.x += bt->radius * (i&1 ? 0.5f : -0.5f);
      // newCenter.y += bt->radius * (i&1 ? 0.5f : -0.5f);

      bt->childs[i] = new Btree_t(bt, newCenter, newRadius);

      cpu_tree_nodes++;
    }
}

/* Inserts a point in the btree creating the appropiate childs
 by setting a maximum number of points*/
void insertPoint2(point_t* point, Btree btree, int maxSize)
{
    if(isLeaf(btree))
    {
      if(btree->points.size() < maxSize)
      {
        btree->points.push_back(point);
      }
      else
      {
        createChilds(btree);
        fillChilds(btree, maxSize);
        int idx = childIdx(point, btree);
        insertPoint2(point, btree->childs[idx], maxSize);

      }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
      int idx = childIdx(point, btree);
      insertPoint2(point, btree->childs[idx], maxSize);
    }
}

/** Calculate the radius in each axis and save the max radius of the bounding box */
Vector2D getRadius(Vector2D min, Vector2D max, float *maxRadius)
{
    Vector2D radii;

    radii.x = (max.x - min.x) / 2.0;
    radii.y = (max.y - min.y) / 2.0;
    // radii.z = (max.z - min.z) / 2.0;

    if(radii.x >= radii.y)
    {
        *maxRadius = radii.x;
    }
    else
    {
        *maxRadius = radii.y;
    }

    return radii;
}


/** Calculate the center of the bounding box */
Vector2D getCenter(Vector2D min, Vector2D radius)
{
    Vector2D center;

    center.x = min.x + radius.x;
    center.y = min.y + radius.y;
    // center.z = min.z + radius.z;

    return center;
}

/* Recursively deletes the tree created by CPU */
void deleteBtree(Btree btree)
{
    if(isLeaf(btree))
    {

      btree->points.clear();
      // btree->childs.clear();

    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 2; i++) {
            // Check
            deleteBtree(btree->childs[i]);
            delete(btree->childs[i]);
        }
        // btree->childs.clear();
    }

    return;
}

int save_areas(std::string file_name, std::string map_name, std::vector<std::pair<Vector2D,float>>& areas)
{
  std::fstream out(file_name, out.out | out.app);

  if(!out.is_open())
    return -1;

  out << map_name << std::endl;

  out.precision(2);
  for(auto& item : areas)
    out << std::fixed << item.first.x << " " << item.first.y << " " << item.second << std::endl;

  out << "\n\n\n";

  out.close();
  if(out.is_open())
    return -1;
  return 0;
}
