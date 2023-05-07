#ifndef TRAVERSE_H

#define TRAVERSE_H

#include <sycl_cuda_ordered/source/sycl_builder.h>

/*Esta estructura se utiliza con puntos ordenados*/

// struct LBVHo
// {
//     const bintree_node*     node_list;
//     const aabb_t*           aabb_list;
//     const point_t*          point_list;
//     uint32_t                numInternalNodes;

//     __device__ __forceinline__ uint32_t  getRoot         (void)           { return 0; }

//     __device__ __forceinline__ bool      isLeaf          (uint32_t node) { return (node >= numInternalNodes); }
//     __device__ __forceinline__ uint32_t  getNumPoints    (uint32_t node) { return node_list[node].numPts; }      
//     __device__ __forceinline__ const aabb_t&  getAABB    (uint32_t node) { return aabb_list[node]; }
//     __device__ __forceinline__ uint32_t  getLeftChild    (uint32_t node) { return node_list[node].left_idx; }
//     __device__ __forceinline__ uint32_t  getRightChild   (uint32_t node) { return node_list[node].right_idx; }

//     __device__ __forceinline__ 
//     point_t getPoint(uint32_t node, uint32_t i) { 
//       int id = node_list[node].object_idx + i;
//       point_t aux = point_list[id];
//       /*aprovecho el cuarto elemento para enviar el índice*/
//       aux.w = real_t(id);
//       return aux; 
//     }

//     __device__ __forceinline__
//     thrust::pair<real_t, uint32_t> getMinPair (uint32_t node) {
//         return thrust::make_pair( point_list[node_list[node].min_idx].z, node_list[node].min_idx );
//     }
// };

struct LBVHo
{
  public:

    LBVHo(const bintree_node* n, const aabb_t* bb, const point_t* p, const uint32_t num):
          node_list(n), aabb_list(bb), point_list(p), numInternalNodes(num) {}

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
    
  private:

    const bintree_node*     node_list;
    const aabb_t*           aabb_list;
    const point_t*          point_list;
    const uint32_t          numInternalNodes;
};


template<typename lbvh_t>
__device__ 
void traverseIterative(lbvh_t& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    uint32_t stack[64];
    uint32_t* stackPtr = stack;
    *stackPtr++ = 0u; // push

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
              }
            }
          }
        }

        __syncthreads();

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
    while (node != 0u);

    // count = numPts;

    return;
}


__global__ void query2D(LBVHo lbvh, aabb_t initBox, double Displace, uint32_t* count, 
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

  if( minNumPoints <= pointsCount){
    count[jdx*nCols + idx] = idmin;
  }

  return;

}



void stage1query2D(SYCL_builder& builder, uint32_t* count, uint32_t Wsize, 
  double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints,
  cl::sycl::queue device_queue){

  // size_t Ncells = nRows*nCols;

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = builder.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.BBox.lower.x + Displace;
  initBox.upper.y = builder.BBox.lower.y + Displace;

  // device_queue.submit([&](cl::sycl::handler &cgh) {

  //   // cgh.parallel_for( cl::sycl::range<1>(nCols*nRows), 
  //   //                   _query(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints));
  //   cgh.parallel_for( cl::sycl::range<2>(nCols,nRows), 
  //                     _query(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints) );

  // }).wait_and_throw();

  // auto mR = sycl::range<1>(nRows*nCols);

  // sycl::buffer<uint32_t, 1> bufMins(count, mR);

  device_queue.submit([&](sycl::handler& h) {

    auto accN = builder.node_list.get_access<sycl::access::mode::read>(h);
    auto accA = builder.aabb_list.get_access<sycl::access::mode::read>(h);
    auto accP = builder.ord_point_cloud.get_access<sycl::access::mode::read>(h);
    // auto accC = bufMins.get_access<sycl::access::mode::write>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dN = reinterpret_cast<node_t*>(ih.get_mem<sycl::backend::cuda>(accN));
      auto dA = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accA));
      auto dP = reinterpret_cast<point_t*>(ih.get_mem<sycl::backend::cuda>(accP));
      // auto dMins = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accC));

      LBVHo lbvh( dN,
                  dA,
                  dP,
                  builder.numInternalNodes );

      // int threadsPerBlock = 64;
      // int numberOfBlocks = int((Ncells-1)/threadsPerBlock) + 1;

      // query<<<numberOfBlocks, threadsPerBlock>>>(lbvh, initBox, Displace, dMins, nCols, Ncells, minNumPoints);

      dim3 tblocks(8,8,1);
      dim3 grid(int((nCols-1)/tblocks.x)+1, int((nRows-1)/tblocks.y)+1, 1);

      query2D<<<grid, tblocks>>>(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints);

    });
  });

  // cudaError_t lastError = cudaGetLastError();
  // if(lastError != cudaSuccess) printf("Error TRAVERSE: %s\n", cudaGetErrorString(lastError));
  // cudaError_t syncError = cudaDeviceSynchronize();
  // if(syncError != cudaSuccess) printf("Sync error TRAVERSE: %s\n", cudaGetErrorString(syncError));

  // auto hMins = bufMins.get_access<sycl::access::mode::read>();

  return;
}

// void stage1query2D(SYCL_builder& builder, uint32_t* count, uint32_t Wsize, 
//   double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints){

//   // size_t Ncells = nRows*nCols;

//   double Displace = round2d(Wsize*(1-Overlap));

//   LBVHo lbvh( builder.node_list,
//               builder.aabb_list,
//               builder.ord_point_cloud,
//               builder.numInternalNodes );

//   aabb_t initBox;
//   initBox.lower.x = builder.BBox.lower.x - Wsize + Displace;
//   initBox.lower.y = builder.BBox.lower.y - Wsize + Displace;
//   initBox.upper.x = builder.BBox.lower.x + Displace;
//   initBox.upper.y = builder.BBox.lower.y + Displace;

//   dim3 tblocks(8,8,1);
//   dim3 grid(int((nCols-1)/tblocks.x)+1, int((nRows-1)/tblocks.y)+1, 1);

// #ifdef DEBUG
//   printf("%d,%d\n", grid.x, grid.y);
// #endif

//   query2D<<<grid, tblocks>>>(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints);

//   return;
// }


#endif // TRAVERSE_H