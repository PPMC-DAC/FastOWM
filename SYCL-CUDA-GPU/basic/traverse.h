#ifndef TRAVERSE_H

#define TRAVERSE_H

#include <bintree/cuda/bintree_builder.h>
#include <ordered/cuda/ordered_builder.h>
#include <octree/cuda/octree_builder.h>

#include <tbb/blocked_range2d.h>

struct LBVH
{
    const bintree_node*   node_list;
    const aabb_t*         aabb_list;
    const uint32_t*       objectOrder;
    const point_t*        point_list;
    uint32_t              numInternalNodes;

    __device__ __forceinline__ uint32_t  getRoot         (void)           { return 0; }

    __device__ __forceinline__ bool      isLeaf          (uint32_t node) { return (node >= numInternalNodes); }
    __device__ __forceinline__ uint32_t  getNumPoints    (uint32_t node) { return node_list[node].numPts; }      
    __device__ __forceinline__ const aabb_t&   getAABB   (uint32_t node) { return aabb_list[node]; }
    __device__ __forceinline__ uint32_t  getLeftChild    (uint32_t node) { return node_list[node].left_idx; }
    __device__ __forceinline__ uint32_t  getRightChild   (uint32_t node) { return node_list[node].right_idx; }
    __device__ __forceinline__ uint32_t  getMinIdx       (uint32_t node) { return node_list[node].min_idx; }

    __device__ __forceinline__ 
    point_t getPoint(uint32_t node, uint32_t i) { 
      point_t aux = point_list[objectOrder[node_list[node].object_idx + i]];
      /*aprovecho el cuarto elemento para enviar el índice*/
      aux.w = real_t(objectOrder[node_list[node].object_idx + i]);
      return aux; 
    }

    __device__ __forceinline__
    thrust::pair<real_t, uint32_t> getMinPair (uint32_t node) {
      return thrust::make_pair( point_list[node_list[node].min_idx].z, node_list[node].min_idx );
    }
};

/*Esta es la estructura para cuando tengo morton64*/
struct LBVHm
{
    const bintree_node*     node_list;
    const aabb_t*           aabb_list;
    const uint64_t*         morton;
    uint32_t                numNodes; /*internal nodes*/

    whole_t whole;
    point_t diffBox;

    __device__ __forceinline__ uint32_t  getRoot         (void)           { return 0; }

    __device__ __forceinline__ bool      isLeaf          (uint32_t node) { return (node >= numNodes); }
    __device__ __forceinline__ uint32_t  getNumPoints    (uint32_t node) { return node_list[node].numPts; }      
    __device__ __forceinline__ const aabb_t&   getAABB    (uint32_t node) { return aabb_list[node]; }
    __device__ __forceinline__ uint32_t  getLeftChild    (uint32_t node) { return node_list[node].left_idx; }
    __device__ __forceinline__ uint32_t  getRightChild   (uint32_t node) { return node_list[node].right_idx; }

    __device__ __forceinline__ 
    point_t getPoint(uint32_t node, uint32_t i) { 
        /*tengo que decodificar el punto porque nos índices apuntan al vector
        de códigos Morton de 64 bits*/
        point_t aux = point3D(morton[node_list[node].object_idx + i]);
        aux.x = aux.x * diffBox.x + whole.lower.x;
        aux.y = aux.y * diffBox.y + whole.lower.y;
        aux.z = aux.z * diffBox.z + whole.lower.z;

        /*aprovecho el cuarto elemento para enviar el índice*/
        aux.w = real_t(node_list[node].object_idx + i);
        return aux; 
    }


    __device__ __forceinline__
    thrust::pair<real_t, uint32_t> getMinPair (uint32_t node) {

        point_t aux = point3D(morton[node_list[node].min_idx]);

        return thrust::make_pair( aux.z * diffBox.z + whole.lower.z, node_list[node].min_idx );
    }
};

/*Esta estructura se utiliza con puntos ordenados*/
struct LBVHo
{
    const bintree_node*     node_list;
    const aabb_t*           aabb_list;
    const point_t*          point_list;
    uint32_t                numInternalNodes;

    __device__ __forceinline__ uint32_t  getRoot         (void)           { return 0; }

    __device__ __forceinline__ bool      isLeaf          (uint32_t node) { return (node >= numInternalNodes); }
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
    thrust::pair<real_t, uint32_t> getMinPair (uint32_t node) {
        return thrust::make_pair( point_list[node_list[node].min_idx].z, node_list[node].min_idx );
    }
};

/* Esta es la estructura para el octree */
struct LBVHoct
{
    const octree_node*    node_list;
    const aabb_t*         aabb_list;
    const uint32_t*       objectOrder;
    const point_t*        point_list;

    __device__ __forceinline__ constexpr uint32_t  getRoot         (void) const           { return 0u; }

    __device__ __forceinline__ constexpr bool      isLeaf          (const uint32_t node) const { return node_list[node].bitmask == 0u; }

    __device__ __forceinline__ constexpr const aabb_t&   getAABB   (const uint32_t node) const { return aabb_list[node]; }

    __device__ __forceinline__ constexpr uint32_t  getChild        (const uint32_t node, const uint32_t i) const { return node_list[node].get_child(i); }
        
    __device__ __forceinline__ 
    const thrust::tuple<const point_t*, uint32_t, uint32_t> getPointPack(const uint32_t node) const { 
      return thrust::make_tuple( &point_list[ node_list[node].object_idx ], node_list[node].object_idx, node_list[node].numPts );
    }
    
    __device__ __forceinline__ 
    const thrust::tuple<real_t, uint32_t, uint32_t> getMinPack(const uint32_t node) const { 
      return thrust::make_tuple( point_list[node_list[node].min_idx].z, node_list[node].min_idx, node_list[node].numPts );
    }
};


/* Esta es la estructura para el octree */
struct LBVHoctCPU
{

  LBVHoctCPU( const octree_node* n, 
              const aabb_t* bb,
              const point_t* p):
              node_list(n), aabb_list(bb), point_list(p) {}

  __host__ __forceinline__ uint32_t  getRoot         (void) const           { return 0u; }
  __host__ __forceinline__ bool      isLeaf          (const uint32_t node) const { return node_list[node].bitmask == 0u; }

  __host__ __forceinline__ const aabb_t&   getAABB   (const uint32_t node) const { return aabb_list[node]; }
  __host__ __forceinline__ uint32_t  getChild        (const uint32_t node, const uint32_t i) const { return node_list[node].get_child(i); }

  
  __host__ __forceinline__ 
  const std::tuple<const point_t*, uint32_t, uint32_t> getPointPack(const uint32_t node) const { 
    return std::make_tuple( &point_list[ node_list[node].object_idx ], node_list[node].object_idx, node_list[node].numPts );
  }
  
  __host__ __forceinline__ 
  const std::tuple<real_t, uint32_t, uint32_t> getMinPack(const uint32_t node) const { 
    return std::make_tuple( point_list[node_list[node].min_idx].z, node_list[node].min_idx, node_list[node].numPts );
  }

  const octree_node*   node_list;
  const aabb_t*   aabb_list;
  const point_t*   point_list;
};



template<typename lbvh_t>
__device__ 
void traverseIterative(lbvh_t& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
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
          auto min = lbvh.getMinPair(childL);
          if(min.first < zmin){
            zmin = min.first;
            idmin = min.second;
          }
        }
        if(isAllOverlaped(queryAABB, lbvh.getAABB(childR))){
          overlapR = false;
          numPts += lbvh.getNumPoints(childR);
          auto min = lbvh.getMinPair(childR);
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
              if(p.z < zmin){
                zmin = p.z;
                idmin = uint32_t(p.w);
                // if(fabs(zmin - lbvh.point_list[idmin].z) > TOL) printf("ERROR\n");
              }
            }
          }
        }
        if (overlapR && lbvh.isLeaf(childR)){
          const uint32_t l = lbvh.getNumPoints(childR);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childR, i);
            if(insideBox(queryAABB, p)){
              numPts++;
              if(p.z < zmin){
                zmin = p.z;
                idmin = uint32_t(p.w);
                // if(fabs(zmin - lbvh.point_list[idmin].z) > TOL) printf("ERROR\n");
              }
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

template<typename lbvh_t>
__device__ 
void countIterative(lbvh_t& lbvh, aabb_t& queryAABB, uint32_t& numPts)
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
        }
        if(isAllOverlaped(queryAABB, lbvh.getAABB(childR))){
          overlapR = false;
          numPts += lbvh.getNumPoints(childR);
        }

        // Query overlaps a leaf node => report collision.

        if (overlapL && lbvh.isLeaf(childL)){
          const uint32_t l = lbvh.getNumPoints(childL);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childL, i);
            if(insideBox(queryAABB, p)){
              numPts++;
            }
          }
        }
        if (overlapR && lbvh.isLeaf(childR)){
          const uint32_t l = lbvh.getNumPoints(childR);
          for(int i=0; i<l; i++){
            point_t p = lbvh.getPoint(childR, i);
            if(insideBox(queryAABB, p)){
              numPts++;
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

//Now using this with Octree_builder and LBVHoct data types
__device__ 
void traverseIterative(LBVHoct& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    uint32_t stack[64];
    uint32_t* stackPtr = stack;
    *stackPtr++ = NULL; // push

    real_t zmin = inf;

    // uint32_t pointsCount = 0;

    // Traverse nodes starting from the root.

    uint32_t node = lbvh.getRoot();
    do
    {
        // Check each child node for overlap.

        // uint32_t childL = lbvh.getLeftChild(node);
        // uint32_t childR = lbvh.getRightChild(node);
        // bool overlapL = (boxOverlap2D(queryAABB, lbvh.getAABB(childL)));
        // bool overlapR = (boxOverlap2D(queryAABB, lbvh.getAABB(childR)));
        uint32_t children[8];
        uint8_t overlap = 0u;

        for(int i=0; i<8; i++){
          children[i] = lbvh.getChild(node,i);
          // overlap |= (children[i] != uint32_t(-1))? (boxOverlap2D(queryAABB, lbvh.getAABB(children[i])) << i) : 0u;
        }

        // Internal node totally overlaped
        // #pragma unroll 8
        for(int i=0; i<8; i++){

          uint32_t child = children[i];

          if(child != uint32_t(-1)){
#ifndef NOMEMO
            if(isAllOverlaped(queryAABB, lbvh.getAABB(child))){
              // overlap = false;

              // numPts += lbvh.getNumPoints(child);
              // auto min = lbvh.getMinPair(child);

              const auto pack = lbvh.getMinPack(child);
              real_t min = thrust::get<0>(pack);
              numPts += thrust::get<2>(pack);

              if(min < zmin){
                zmin = min;
                idmin = thrust::get<1>(pack);
              }
            }
            else
            {
#endif
              overlap |= boxOverlap2D(queryAABB, lbvh.getAABB(child)) << i;
#ifndef NOMEMO
            }
#endif
          }
        }
        // __syncthreads();
        // Query overlaps a leaf node => report collision.
        #pragma unroll 8
        for(int i=0; i<8; i++){

          uint32_t child = children[i];

          if ( (overlap & (1u << i)) && lbvh.isLeaf(child) ){

            // const uint32_t l = lbvh.getNumPoints(child);
            // const point_t* p = lbvh.getPointList(child);
            // const uint32_t pointID = lbvh.getObjID(child);

            const auto pack = lbvh.getPointPack(child);
            const point_t* p = thrust::get<0>(pack);
            // const uint32_t pointID = thrust::get<1>(pack);
            const uint32_t l = thrust::get<2>(pack);

            for(int i=0; i<l; i++){
              // point_t p = lbvh.getPoint(child, i);
              if(insideBox(queryAABB, p[i])){
                numPts++;
                if(p[i].z < zmin){
                  zmin = p[i].z;
                  idmin = thrust::get<1>(pack) + i;
                }
              }
            }
          }
        }
        __syncthreads();

        // Query overlaps an internal node => traverse.

        uint8_t traverse = 0u;
        for(int i=0; i<8; i++){
          traverse |= ( (overlap & (1u << i)) && !lbvh.isLeaf(children[i]) ) << i;
        }

        for(int i=0; i<8; i++)
        {
          if(traverse & (1u << i))
            *stackPtr++ = children[i]; // push
        }
        node = *--stackPtr; // pop

        // __syncthreads();
    }
    while (node != NULL);

    // numPts = pointsCount;

    return;
}

__host__
void traverseIterativeCPU(const LBVHoctCPU& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin, uint32_t isTmp)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    uint32_t stack[64];
    uint32_t* stackPtr = stack;
    *stackPtr++ = 0u; // push

    real_t zmin = inf;

    uint32_t node = lbvh.getRoot();
    
    do
    {

      for(int i=0; i<8; i++){

        uint32_t child = lbvh.getChild(node,i);
        // int overlap = 0;

        if(child != uint32_t(-1)){
          
          if( isAllOverlaped(queryAABB, lbvh.getAABB(child)) ){ // completely overlaped

            const auto pack = lbvh.getMinPack(child);
            real_t min = std::get<0>(pack);
            numPts += std::get<2>(pack);

            if(min < zmin){
              zmin = min;
              idmin = std::get<1>(pack);
            }
          }
          else
          {
            // overlap = boxOverlap2D(queryAABB, lbvh.getAABB(child));

            if( boxOverlap2D(queryAABB, lbvh.getAABB(child)) ) // overlaped
            {

              if ( lbvh.isLeaf(child) ) { // leaf

                const auto pack = lbvh.getPointPack(child);
                const point_t* p = std::get<0>(pack);
                const uint32_t l = std::get<2>(pack);

                for(int i=0; i<l; i++){
                  if(insideBox(queryAABB, p[i])){
                    numPts++;
                    if(p[i].z < zmin){
                      zmin = p[i].z;
                      idmin = std::get<1>(pack) + i;
                    }
                  }
                }

              } else { // no leaf

                *stackPtr++ = child; //push

              }

            }

          }

        } // child

      } // for

      node = *--stackPtr; // pop

    }
    while (node != 0u);

    // if(isTmp && idmin == 0xFFFFFFFFu){
    //   printf("Minimo: %u; PROBLEMMMMM: %u\n", idmin, 0xFFFFFFFFu);
    // } 

    return;
}


__global__ void query(const bintree_node* n, const aabb_t* bb, const uint64_t* m,
    aabb_t initBox, whole_t BBox, point_t diffBox, double Displace, uint32_t* count, 
    uint32_t numNodes, uint32_t numCells, uint32_t minNumPoints, uint32_t nCols)
{
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx >= numCells)
      return;

  LBVHm lbvh;
  lbvh.node_list       = n;
  lbvh.aabb_list       = bb;
  lbvh.morton          = m;
  lbvh.numNodes        = numNodes;

  lbvh.whole = BBox;
  lbvh.diffBox = diffBox;

  aabb_t cellBox;
  cellBox.upper.x = initBox.upper.x + (idx%nCols)*Displace;
  cellBox.lower.x = initBox.lower.x + (idx%nCols)*Displace;
  cellBox.upper.y = initBox.upper.y + (int)(idx/nCols)*Displace;
  cellBox.lower.y = initBox.lower.y + (int)(idx/nCols)*Displace;

  uint32_t pointsCount = 0;
  uint32_t idmin = 0xFFFFFFFFu;

  traverseIterative(lbvh, cellBox, pointsCount, idmin);

  if( minNumPoints <= pointsCount){
    count[idx] = idmin;
  }

  return;

}

void stage1query(bintree_node* n, aabb_t* bb, uint64_t* m, uint32_t* count,
  uint32_t Wsize, double Overlap, uint32_t nCols, uint32_t nRows,
  uint32_t minNumPoints, whole_t BBox, point_t diffBox, uint32_t numNodes){

    size_t Ncells = nRows*nCols;

    double Displace = round2d(Wsize*(1-Overlap));

    aabb_t initBox;
    initBox.lower.x = BBox.lower.x - Wsize + Displace;
    initBox.lower.y = BBox.lower.y - Wsize + Displace;
    initBox.upper.x = BBox.lower.x + Displace;
    initBox.upper.y = BBox.lower.y + Displace;

    int threadsPerBlock = 64;
    int numberOfBlocks = int(Ncells/threadsPerBlock - 1.0) + 1;

    query<<<numberOfBlocks, threadsPerBlock>>>(n, bb, m, initBox, BBox, diffBox, Displace, count, numNodes, Ncells, minNumPoints, nCols);

    cudaError_t lastError = cudaGetLastError();
    if(lastError != cudaSuccess) printf("Error KERNEL: %s\n", cudaGetErrorString(lastError));
    cudaDeviceSynchronize();


  return;
}

template<typename lbvh_t>
__global__ void query(lbvh_t lbvh, aabb_t initBox, double Displace, uint32_t* count, 
    uint32_t nCols, uint32_t numCells, uint32_t minNumPoints)
{
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx >= numCells)
      return;

  aabb_t cellBox;
  cellBox.upper.x = initBox.upper.x + (idx%nCols)*Displace;
  cellBox.lower.x = initBox.lower.x + (idx%nCols)*Displace;
  cellBox.upper.y = initBox.upper.y + (int)(idx/nCols)*Displace;
  cellBox.lower.y = initBox.lower.y + (int)(idx/nCols)*Displace;

#ifdef DEBUG
  if(idx == numCells-1){
    printf("%d,%d\n", (idx%nCols), (int)(idx/nCols));
    printf("%lf, %lf  %lf, %lf\n", cellBox.upper.x, cellBox.upper.y, cellBox.lower.x, cellBox.lower.y);
  }
#endif

  uint32_t pointsCount = 0;
  uint32_t idmin = 0xFFFFFFFFu;

  traverseIterative(lbvh, cellBox, pointsCount, idmin);

  if( minNumPoints <= pointsCount){
    count[idx] = idmin;
  }

  return;

}

void stage1query(Ordered_builder& builder, uint32_t* count, uint32_t Wsize, 
  double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints){

  size_t Ncells = nRows*nCols;

  double Displace = round2d(Wsize*(1-Overlap));

  LBVHo lbvh;
  lbvh.node_list        = builder.node_list;
  lbvh.aabb_list        = builder.aabb_list;
  lbvh.point_list       = builder.ord_point_cloud;
  lbvh.numInternalNodes = builder.numInternalNodes;

  aabb_t initBox;
  initBox.lower.x = builder.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.BBox.lower.x + Displace;
  initBox.upper.y = builder.BBox.lower.y + Displace;

  int threadsPerBlock = 64;
  int numberOfBlocks = int((Ncells-1)/threadsPerBlock) + 1;

  query<<<numberOfBlocks, threadsPerBlock>>>(lbvh, initBox, Displace, count, nCols, Ncells, minNumPoints);

  cudaError_t lastError = cudaGetLastError();
  if(lastError != cudaSuccess) printf("Error KERNEL: %s\n", cudaGetErrorString(lastError));
  cudaDeviceSynchronize();


  return;
}

template<typename lbvh_t>
__global__ void query2D(lbvh_t lbvh, aabb_t initBox, double Displace, uint32_t* count, 
    uint32_t nCols, uint32_t nRows, uint32_t minNumPoints)
{
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  int jdx = threadIdx.y + blockDim.y * blockIdx.y;
  if (idx >= nCols || jdx >= nRows)
      return;
  
  aabb_t cellBox;
  cellBox.upper.x = idx*Displace + initBox.upper.x;
  cellBox.lower.x = idx*Displace + initBox.lower.x;
  cellBox.upper.y = jdx*Displace + initBox.upper.y;
  cellBox.lower.y = jdx*Displace + initBox.lower.y;

#ifdef DEBUG
  if(idx == nCols-1 && jdx == nRows-1){
    printf("%d,%d\n", idx, jdx);
    printf("%lf, %lf  %lf, %lf\n", cellBox.upper.x, cellBox.upper.y, cellBox.lower.x, cellBox.lower.y);
  }
#endif

  uint32_t pointsCount = 0;
  uint32_t idmin = 0xFFFFFFFFu;

  traverseIterative(lbvh, cellBox, pointsCount, idmin);

  // if( minNumPoints <= pointsCount){
  //   count[jdx*nCols + idx] = idmin;
  // }

  count[jdx*nCols + idx] = ( minNumPoints <= pointsCount )? idmin : 0xFFFFFFFFu;

  return;

}

template<typename lbvh_t>
__global__ void queryStage3(lbvh_t lbvh, lbvh_t lbvh_grid, aabb_t initBox, uint32_t Displace, uint32_t* count, 
    uint32_t nCols, uint32_t nRows, uint32_t minNumPoints)
{
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  int jdx = threadIdx.y + blockDim.y * blockIdx.y;
  if (idx >= nCols || jdx >= nRows)
      return;
  
  aabb_t cellBox;
  cellBox.upper.x = idx*Displace + initBox.upper.x;
  cellBox.lower.x = idx*Displace + initBox.lower.x;
  cellBox.upper.y = jdx*Displace + initBox.upper.y;
  cellBox.lower.y = jdx*Displace + initBox.lower.y;

#ifdef DEBUG
  if(idx == nCols-1 && jdx == nRows-1){
    printf("GRID: %d,%d\n", idx, jdx);
    printf("%lf, %lf  %lf, %lf\n", cellBox.upper.x, cellBox.upper.y, cellBox.lower.x, cellBox.lower.y);
  }
#endif

  uint32_t pointsCount = 0;
  uint32_t idmin = 0xFFFFFFFFu;

  countIterative(lbvh_grid, cellBox, pointsCount);

  if(pointsCount == 0)
  {
    traverseIterative(lbvh, cellBox, pointsCount, idmin);

    count[jdx*nCols + idx] = ( pointsCount > 0 )? idmin : 0xFFFFFFFFu;

  }else{
    count[jdx*nCols + idx] = 0xFFFFFFFFu;
  }

  return;

}

void stage1query2D(Ordered_builder& builder, uint32_t* count, uint32_t Wsize, 
  double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints){

  // size_t Ncells = nRows*nCols;

  double Displace = round2d(Wsize*(1-Overlap));

  LBVHo lbvh;
  lbvh.node_list        = builder.node_list;
  lbvh.aabb_list        = builder.aabb_list;
  lbvh.point_list       = builder.ord_point_cloud;
  lbvh.numInternalNodes = builder.numInternalNodes;

  aabb_t initBox;
  initBox.lower.x = builder.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.BBox.lower.x + Displace;
  initBox.upper.y = builder.BBox.lower.y + Displace;

  // int threadsPerBlock = 64;
  // int numberOfBlocks = int(Ncells/threadsPerBlock - 1.0) + 1;
  // query<<<numberOfBlocks, threadsPerBlock>>>(lbvh, initBox, Displace, count, nCols, Ncells, minNumPoints);

  dim3 tblocks(8,8,1);
  dim3 grid(int((nCols-1)/tblocks.x)+1, int((nRows-1)/tblocks.y)+1, 1);
#ifdef DEBUG
  printf("%d,%d\n", grid.x, grid.y);
#endif

  query2D<<<grid, tblocks>>>(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints);

  cudaError_t lastError = cudaGetLastError();
  if(lastError != cudaSuccess) printf("Error KERNEL S1: %s\n", cudaGetErrorString(lastError));
  lastError =cudaDeviceSynchronize();
  if(lastError != cudaSuccess) printf("Error KERNEL S1 SYNC: %s\n", cudaGetErrorString(lastError));

  return;
}

void stage3query2D(Ordered_builder& builder, Ordered_builder& builder_grid, uint32_t* count, uint32_t Bsize,
 uint32_t nCols, uint32_t nRows, uint32_t minNumPoints){

  // size_t Ncells = nRows*nCols;

  LBVHo lbvh;
  lbvh.node_list        = builder.node_list;
  lbvh.aabb_list        = builder.aabb_list;
  lbvh.point_list       = builder.ord_point_cloud;
  lbvh.numInternalNodes = builder.numInternalNodes;

  LBVHo lbvh_grid;
  lbvh_grid.node_list        = builder_grid.node_list;
  lbvh_grid.aabb_list        = builder_grid.aabb_list;
  lbvh_grid.point_list       = builder_grid.ord_point_cloud;
  lbvh_grid.numInternalNodes = builder_grid.numInternalNodes;

  aabb_t initBox;
  initBox.lower.x = builder.BBox.lower.x;
  initBox.lower.y = builder.BBox.lower.y;
  initBox.upper.x = builder.BBox.lower.x + Bsize;
  initBox.upper.y = builder.BBox.lower.y + Bsize;

  // int threadsPerBlock = 64;
  // int numberOfBlocks = int(Ncells/threadsPerBlock - 1.0) + 1;
  // query<<<numberOfBlocks, threadsPerBlock>>>(lbvh, initBox, Displace, count, nCols, Ncells, minNumPoints);

  dim3 tblocks(8,8,1);
  dim3 grid(int((nCols-1)/tblocks.x)+1, int((nRows-1)/tblocks.y)+1, 1);
#ifdef DEBUG
  printf("%d,%d\n", grid.x, grid.y);
#endif

  queryStage3<<<grid, tblocks>>>(lbvh, lbvh_grid, initBox, Bsize, count, nCols, nRows, minNumPoints);

  cudaError_t lastError = cudaGetLastError();
  if(lastError != cudaSuccess) printf("Error KERNEL S3: %s\n", cudaGetErrorString(lastError));
  lastError =cudaDeviceSynchronize();
  if(lastError != cudaSuccess) printf("Error KERNEL S3 SYNC: %s\n", cudaGetErrorString(lastError));
  
  return;
}

//Now using this with Octree_builder data type
void stage1query2D(Octree_builder& builder, uint32_t* count, uint32_t Wsize, 
  double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints){

  double Displace = round2d(Wsize*(1-Overlap));

  // LBVHoct lbvh;
  // lbvh.node_list        = thrust::raw_pointer_cast( &builder.m_octree.front() );
  // lbvh.aabb_list        = thrust::raw_pointer_cast( &builder.m_aabb.front() );
  // lbvh.point_list       = builder.bintree.ord_point_cloud;
  LBVHoct lbvh;
  lbvh.node_list        = builder.b_octree;
  lbvh.aabb_list        = builder.b_aabb;
  lbvh.point_list       = builder.bintree.ord_point_cloud;

  aabb_t initBox;
  initBox.lower.x = builder.bintree.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.bintree.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.bintree.BBox.lower.x + Displace;
  initBox.upper.y = builder.bintree.BBox.lower.y + Displace;

  dim3 tblocks(8,8,1);
  dim3 grid(int((nCols-1)/tblocks.x)+1, int((nRows-1)/tblocks.y)+1, 1);

  query2D<<<grid, tblocks>>>(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints);

  // cudaError_t lastError = cudaGetLastError();
  // if(lastError != cudaSuccess) printf("Error KERNEL: %s\n", cudaGetErrorString(lastError));
  // cudaDeviceSynchronize();


  return;
}

/*En este caso pongo un template para que no haya problemas cuando empiezo
desde cero con una nueva clase. Como utilizo como plantilla el caso de 
"bintree", pongo el template en el método para el recorrido con "Bintree_builder"*/
template<typename node_builder>
void stage1query2D(node_builder& builder, uint32_t* count, uint32_t Wsize, 
  double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints){

  double Displace = round2d(Wsize*(1-Overlap));

  LBVH lbvh;
  lbvh.node_list        = builder.node_list;
  lbvh.aabb_list        = builder.aabb_list;
  lbvh.objectOrder      = builder.indices_out;
  lbvh.point_list       = builder.point_cloud;
  lbvh.numInternalNodes = builder.numInternalNodes;

  aabb_t initBox;
  initBox.lower.x = builder.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.BBox.lower.x + Displace;
  initBox.upper.y = builder.BBox.lower.y + Displace;

  dim3 tblocks(8,8,1);
  dim3 grid(int((nCols-1)/tblocks.x)+1, int((nRows-1)/tblocks.y)+1, 1);

  query2D<<<grid, tblocks>>>(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints);

  cudaError_t lastError = cudaGetLastError();
  if(lastError != cudaSuccess) printf("Error KERNEL: %s\n", cudaGetErrorString(lastError));
  cudaDeviceSynchronize();


  return;
}


class _query{
  
  public:
    _query(const LBVHoctCPU _lbvh, const aabb_t _initBox, const double _Displace, const double _Overlap,
            uint32_t* _count, const uint32_t _nCols, const uint32_t _minNumPoints) : 
          lbvh(_lbvh), initBox(_initBox), Displace(_Displace), Overlap(_Overlap),
          count(_count), nCols(_nCols), minNumPoints(_minNumPoints) {}
    
    void operator()(tbb::blocked_range2d<int,int> r) const
    {
      // int idx = static_cast<int>(index[0]);
      // int jdx = static_cast<int>(index[1]);
      int je = r.rows().end();
      int ie = r.cols().end();
      
      aabb_t cellBox;

      for(int jj = r.rows().begin(); jj < je; jj++ ) {

        uint32_t pointsCount;
        uint32_t idmin = 0xFFFFFFFFu;

        // cellCenter.y = initY + jj*Displace;
        cellBox.upper.y = jj*Displace + initBox.upper.y;
        cellBox.lower.y = jj*Displace + initBox.lower.y;

        for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

          // cellCenter.x = initX + ii*Displace;
          cellBox.upper.x = ii*Displace + initBox.upper.x;
          cellBox.lower.x = ii*Displace + initBox.lower.x;

          if(idmin != 0xFFFFFFFFu && insideBox(cellBox, lbvh.point_list[idmin]))
          {

            // Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};
            aabb_t oCell;
            oCell.upper.x = cellBox.upper.x;
            oCell.lower.x = cellBox.upper.x - Displace;
            oCell.upper.y = cellBox.upper.y;
            oCell.lower.y = cellBox.lower.y;

            uint32_t old_cellPoints = pointsCount;

            pointsCount = 0u;
            uint32_t tmpId = 0xFFFFFFFFu;

            traverseIterativeCPU(lbvh, oCell, pointsCount, tmpId, 1u);

            // We're assuming the points were equidistant throughout the cell, which isn't always true.

            /*En este punto, si queremos ser estrictos, en vez de hacer esta suposición podemos
            lanzar una búsqueda "countMin" en la zona en la que conocemos el mínimo, pudiendo lanzar
            las dos búsquedas diferentes en paralelo con un "parallel_invoke" */
            pointsCount += (uint32_t)(old_cellPoints * Overlap);

            /*Si he hecho traverseIterativeCPU y tmpId = 0xFFFFFFFFu, quiere decir que la BB con
            la que analizo no solapa ninguna BB de un nodo, lo que puede suceder muy facilmente
            porque la ventana de análisis tiene un margen de tamaño (Wsize-Displace) por cada
            lateral de la nube de puntos*/
            if(tmpId != 0xFFFFFFFFu && lbvh.point_list[tmpId].z < lbvh.point_list[idmin].z)
            {
              idmin = tmpId;
            }

          } else {

            pointsCount = 0u;
            idmin = 0xFFFFFFFFu;

            traverseIterativeCPU(lbvh, cellBox, pointsCount, idmin, 0u);

          }

          count[jj*nCols + ii] = ( minNumPoints <= pointsCount )? idmin : 0u;

        }

      }


      return;
    }
  
  private:
    const LBVHoctCPU  lbvh;
    const aabb_t  initBox;
    const double Overlap;
    const double Displace;
    uint32_t* count;
    const uint32_t nCols;
    // const uint32_t nRows;
    const uint32_t minNumPoints;

};

void stage1query2DCPU(const Octree_builder& builder, uint32_t* count, const uint32_t Wsize, 
  const double Overlap, const uint32_t wCols, const uint32_t nCols, const uint32_t nRows, const uint32_t minNumPoints){

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = builder.bintree.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.bintree.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.bintree.BBox.lower.x + Displace;
  initBox.upper.y = builder.bintree.BBox.lower.y + Displace;


  LBVHoctCPU lbvh(
    builder.b_octree,
    builder.b_aabb,
    builder.bintree.ord_point_cloud
  );

  // const size_t dimChunk = 8;

  tbb::parallel_for( tbb::blocked_range2d<int,int>{ 0, static_cast<int>(nRows),
                                                    static_cast<int>(wCols) ,static_cast<int>(nCols)},
                    _query(lbvh, initBox, Displace, Overlap, count, nCols, minNumPoints) );

  return;
}

#endif // TRAVERSE_H