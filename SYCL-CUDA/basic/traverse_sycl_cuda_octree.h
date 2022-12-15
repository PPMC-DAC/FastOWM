#ifndef TRAVERSE_H

#define TRAVERSE_H

// #include <bintree/cuda/bintree_builder.h>
// #include <ordered/cuda/ordered_builder.h>
// #include <octree/cuda/octree_builder.h>
#include <sycl_cuda_octree/cuda/octree_builder.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>


/* Esta es la estructura para el octree */
struct LBVHoct
{
    const octree_node*    node_list;
    const aabb_t*         aabb_list;
    // const uint32_t*       objectOrder;
    const point_t*        point_list;

    __device__ __forceinline__ uint32_t  getRoot         (void)           { return 0u; }
    __device__ __forceinline__ bool      isLeaf          (const uint32_t node) { return node_list[node].bitmask == 0u; }

    __device__ __forceinline__ const aabb_t&   getAABB   (const uint32_t node) { return aabb_list[node]; }
    __device__ __forceinline__ uint32_t  getChild        (const uint32_t node, const uint32_t i) { return node_list[node].get_child(i); }

    
    __device__ __forceinline__ 
    const thrust::tuple<const point_t*, uint32_t, uint32_t> getPointPack(const uint32_t node) { 
      return thrust::make_tuple( &point_list[ node_list[node].object_idx ], node_list[node].object_idx, node_list[node].numPts );
    }
    
    __device__ __forceinline__ 
    const thrust::tuple<real_t, uint32_t, uint32_t> getMinPack(const uint32_t node) { 
      return thrust::make_tuple( point_list[node_list[node].min_idx].z, node_list[node].min_idx, node_list[node].numPts );
    }
};

/* Esta es la estructura para el octree */
template<typename accN_t, typename accA_t, typename accP_t>
struct LBVHoctCPU
{

  LBVHoctCPU( accN_t n, 
              accA_t bb,
              accP_t p):
        node_list(n), aabb_list(bb), point_list(p) {}

  inline uint32_t  getRoot         (void) const           { return 0u; }
  inline bool      isLeaf          (const uint32_t node) const { return node_list[node].bitmask == 0u; }

  inline const aabb_t&   getAABB   (const uint32_t node) const { return aabb_list[node]; }
  inline uint32_t  getChild        (const uint32_t node, const uint32_t i) const { return node_list[node].get_child(i); }

  
  inline 
  const std::tuple<const point_t*, uint32_t, uint32_t> getPointPack(const uint32_t node) const { 
    return std::make_tuple( &point_list[ node_list[node].object_idx ], node_list[node].object_idx, node_list[node].numPts );
  }
  
  inline 
  const std::tuple<real_t, uint32_t, uint32_t> getMinPack(const uint32_t node) const { 
    return std::make_tuple( point_list[node_list[node].min_idx].z, node_list[node].min_idx, node_list[node].numPts );
  }

  accN_t   node_list;
  accA_t   aabb_list;
  accP_t   point_list;
};


__device__
void traverseIterative(LBVHoct& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    uint32_t stack[64];
    uint32_t* stackPtr = stack;
    *stackPtr++ = 0u; // push

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
              overlap |= boxOverlap2D(queryAABB, lbvh.getAABB(child)) << i;
            }
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
    while (node != 0u);

    // numPts = pointsCount;

    return;
}

template<typename lbvh_t>
void traverseIterativeCPU(lbvh_t& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin, uint32_t isTmp)
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

template<typename lbvh_t>
class _query{
  
  public:
    _query(lbvh_t _lbvh, const aabb_t _initBox, const double _Displace, const double _Overlap,
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
    lbvh_t  lbvh;
    const aabb_t  initBox;
    const double Overlap;
    const double Displace;
    uint32_t* count;
    const uint32_t nCols;
    // const uint32_t nRows;
    const uint32_t minNumPoints;

};



__global__ void query(LBVHoct lbvh, aabb_t initBox, double Displace, uint32_t* count, 
    uint32_t nCols, uint32_t Ncells, uint32_t minNumPoints)
{
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx >= Ncells)
      return;
  
  aabb_t cellBox;
  cellBox.upper.x = initBox.upper.x + (idx%nCols)*Displace;
  cellBox.lower.x = initBox.lower.x + (idx%nCols)*Displace;
  cellBox.upper.y = initBox.upper.y + (int)(idx/nCols)*Displace;
  cellBox.lower.y = initBox.lower.y + (int)(idx/nCols)*Displace;

// #ifdef DEBUG
//   if(idx == nCols-1 && jdx == nRows-1){
//     printf("%d,%d\n", idx, jdx);
//     printf("%lf, %lf  %lf, %lf\n", cellBox.upper.x, cellBox.upper.y, cellBox.lower.x, cellBox.lower.y);
//   }
// #endif

  uint32_t pointsCount = 0;
  uint32_t idmin = 0xFFFFFFFFu;

  traverseIterative(lbvh, cellBox, pointsCount, idmin);

  // if( minNumPoints <= pointsCount){
  //   count[jdx*nCols + idx] = idmin;
  // }

  count[idx] = ( minNumPoints <= pointsCount )? idmin : 0u;

  return;

}

__global__ void query2D(LBVHoct lbvh, aabb_t initBox, double Displace, uint32_t* count, 
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

  count[jdx*nCols + idx] = ( minNumPoints <= pointsCount )? idmin : 0u;

  return;

}



void stage1query(Octree_builder& builder, sycl::buffer<uint32_t, 1> b_count, uint32_t Wsize, 
  double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints,
  cl::sycl::queue device_queue){

  uint32_t Ncells = nCols*nRows;

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = builder.bintree.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.bintree.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.bintree.BBox.lower.x + Displace;
  initBox.upper.y = builder.bintree.BBox.lower.y + Displace;

  device_queue.submit([&](sycl::handler& h) {

    sycl::accessor accN(builder.b_octree, h, sycl::read_only);
    sycl::accessor accA(builder.b_aabb, h, sycl::read_only);
    sycl::accessor accP(builder.bintree.ord_point_cloud, h, sycl::read_only);
    sycl::accessor accM(b_count, h, sycl::write_only);

    h.interop_task([=](sycl::interop_handler ih) {
      
      auto dN = reinterpret_cast<octree_node*>(ih.get_mem<sycl::backend::cuda>(accN));
      auto dA = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accA));
      auto dP = reinterpret_cast<point_t*>(ih.get_mem<sycl::backend::cuda>(accP));
      auto dMins = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accM));

      LBVHoct lbvh;

      lbvh.node_list = dN;
      lbvh.aabb_list = dA;
      lbvh.point_list = dP;

      int threadsPerBlock = 64;
      int numberOfBlocks = int((Ncells-1)/threadsPerBlock) + 1;

      query<<<numberOfBlocks, threadsPerBlock>>>(lbvh, initBox, Displace, dMins, nCols, Ncells, minNumPoints);

      // dim3 tblocks(8,8,1);
      // dim3 grid(int((nCols-1)/tblocks.x)+1, int((nRows-1)/tblocks.y)+1, 1);

      // query2D<<<grid, tblocks>>>(lbvh, initBox, Displace, dMins, nCols, nRows, minNumPoints);

    });
  });

  return;
}


void stage1query2D(Octree_builder& builder, sycl::buffer<uint32_t, 1> b_count, uint32_t Wsize, 
  double Overlap, uint32_t wCols, uint32_t nRows, uint32_t minNumPoints,
  cl::sycl::queue device_queue){

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = builder.bintree.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.bintree.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.bintree.BBox.lower.x + Displace;
  initBox.upper.y = builder.bintree.BBox.lower.y + Displace;

  device_queue.submit([&](sycl::handler& h) {

    sycl::accessor accN(builder.b_octree, h, sycl::read_only);
    sycl::accessor accA(builder.b_aabb, h, sycl::read_only);
    sycl::accessor accP(builder.bintree.ord_point_cloud, h, sycl::read_only);
    sycl::accessor accM(b_count, h, sycl::write_only);

    h.interop_task([=](sycl::interop_handler ih) {
      
      auto dN = reinterpret_cast<octree_node*>(ih.get_mem<sycl::backend::cuda>(accN));
      auto dA = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accA));
      auto dP = reinterpret_cast<point_t*>(ih.get_mem<sycl::backend::cuda>(accP));
      auto dMins = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accM));

      LBVHoct lbvh;

      lbvh.node_list = dN;
      lbvh.aabb_list = dA;
      lbvh.point_list = dP;

      // int threadsPerBlock = 64;
      // int numberOfBlocks = int((Ncells-1)/threadsPerBlock) + 1;

      // query<<<numberOfBlocks, threadsPerBlock>>>(lbvh, initBox, Displace, dMins, nCols, Ncells, minNumPoints);

      // int wCols = nCols*3/4;

      dim3 tblocks(8,8,1);
      dim3 grid(int((wCols-1)/tblocks.x)+1, int((nRows-1)/tblocks.y)+1, 1);

      query2D<<<grid, tblocks>>>(lbvh, initBox, Displace, dMins, wCols, nRows, minNumPoints);

    });
  });

  return;
}

void stage1query2DCPU(Octree_builder& builder, uint32_t* count, const uint32_t Wsize, 
  const double Overlap, const uint32_t wCols, const uint32_t nCols, const uint32_t nRows, const uint32_t minNumPoints){

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = builder.bintree.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.bintree.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.bintree.BBox.lower.x + Displace;
  initBox.upper.y = builder.bintree.BBox.lower.y + Displace;


  sycl::host_accessor accN(builder.b_octree, sycl::read_only);
  sycl::host_accessor accA(builder.b_aabb, sycl::read_only);
  sycl::host_accessor accP(builder.bintree.ord_point_cloud, sycl::read_only);

  LBVHoctCPU lbvh(
    accN,
    accA,
    accP
  );

  // const size_t dimChunk = 8;

  tbb::parallel_for( tbb::blocked_range2d<int,int>{ 0, static_cast<int>(nRows),
                                                    static_cast<int>(wCols) ,static_cast<int>(nCols)},
                    _query(lbvh, initBox, Displace, Overlap, count, nCols, minNumPoints) );

  return;
}

template<typename lbvh_t>
void stage1query2DCPU(lbvh_t& lbvh, whole_t& whole, uint32_t* count, uint32_t Wsize, 
  double Overlap, uint32_t wCols, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints){

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = whole.lower.x - Wsize + Displace;
  initBox.lower.y = whole.lower.y - Wsize + Displace;
  initBox.upper.x = whole.lower.x + Displace;
  initBox.upper.y = whole.lower.y + Displace;

  // size_t dimChunk = 2;

  tbb::parallel_for( tbb::blocked_range2d<int,int>{ 0, static_cast<int>(nRows),
                                                    static_cast<int>(wCols) ,static_cast<int>(nCols)},
                    _query(lbvh, initBox, Displace, Overlap, count, nCols, minNumPoints) );

  return;
}


#endif // TRAVERSE_H