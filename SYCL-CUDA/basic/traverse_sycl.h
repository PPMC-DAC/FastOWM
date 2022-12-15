#ifndef TRAVERSE_H

#define TRAVERSE_H

// #include <bintree/cuda/bintree_builder.h>
// #include <ordered/cuda/ordered_builder.h>
// #include <octree/cuda/octree_builder.h>
#include <sycl_ordered/source/sycl_builder.h>

#include <tbb/blocked_range2d.h>

/*Esta estructura se utiliza con puntos ordenados*/
struct LBVHo
{
  public:

    LBVHo(const bintree_node* n, const aabb_t* bb, const point_t* p, const uint32_t num):
          node_list(n), aabb_list(bb), point_list(p), numInternalNodes(num) {}

    inline uint32_t  getRoot         (void) const           { return 0; }

    inline bool      isLeaf          (uint32_t node) const { return (node >= numInternalNodes); }
    inline uint32_t  getNumPoints    (uint32_t node) const { return node_list[node].numPts; }      
    inline const aabb_t&  getAABB    (uint32_t node) const { return aabb_list[node]; }
    inline uint32_t  getLeftChild    (uint32_t node) const { return node_list[node].left_idx; }
    inline uint32_t  getRightChild   (uint32_t node) const { return node_list[node].right_idx; }

    inline 
    point_t getPoint(uint32_t node, uint32_t i) const { 
      int id = node_list[node].object_idx + i;
      point_t aux = point_list[id];
      /*aprovecho el cuarto elemento para enviar el índice*/
      aux.w = real_t(id);
      return aux; 
    }

    inline
    std::pair<real_t, uint32_t> getMinPair (uint32_t node) const {
        return std::make_pair( point_list[node_list[node].min_idx].z, node_list[node].min_idx );
    }
    
  private:
    const bintree_node*     node_list;
    const aabb_t*           aabb_list;
    const point_t*          point_list;
    const uint32_t          numInternalNodes;
};

struct LBVHoct
{

  LBVHoct( const octree_node* n, 
              const aabb_t* bb,
              const point_t* p):
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

  const octree_node*   node_list;
  const aabb_t*   aabb_list;
  const point_t*   point_list;
};


inline void traverseIterative(const LBVHo& lbvh, const aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
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

class _query{
  
  public:
    _query(const LBVHo _lbvh, const aabb_t _initBox, const double _Displace, uint32_t* _count, 
          const uint32_t _nCols, const uint32_t _nRows, const uint32_t _minNumPoints) : 
          lbvh(_lbvh), initBox(_initBox), Displace(_Displace), count(_count),
          nCols(_nCols), nRows(_nRows), minNumPoints(_minNumPoints) {}

    void operator()(cl::sycl::id<1> index) const
    {

      int idx = static_cast<int>(index[0]);
      if(idx >= nCols*nRows)
        return;

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
    
    void operator()(cl::sycl::id<2> index) const
    {
      int idx = static_cast<int>(index[0]);
      int jdx = static_cast<int>(index[1]);

      if (idx >= nCols || jdx >= nRows)
          return;
      
      aabb_t cellBox;
      cellBox.upper.x = idx*Displace + initBox.upper.x;
      cellBox.lower.x = idx*Displace + initBox.lower.x;
      cellBox.upper.y = jdx*Displace + initBox.upper.y;
      cellBox.lower.y = jdx*Displace + initBox.lower.y;

      uint32_t pointsCount = 0;
      uint32_t idmin = 0xFFFFFFFFu;

      traverseIterative(lbvh, cellBox, pointsCount, idmin);

      if( minNumPoints <= pointsCount){
        count[jdx*nCols + idx] = idmin;
      }

      return;
    }
  
  private:
    const LBVHo  lbvh;
    const aabb_t  initBox;
    const double Displace;
    uint32_t* count;
    const uint32_t nCols;
    const uint32_t nRows;
    const uint32_t minNumPoints;

};

void stage1query2D(SYCL_builder& builder, uint32_t* count, uint32_t Wsize, 
  double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints,
  cl::sycl::queue device_queue){

  // size_t Ncells = nRows*nCols;

  double Displace = round2d(Wsize*(1-Overlap));

  LBVHo lbvh( builder.node_list,
              builder.aabb_list,
              builder.ord_point_cloud,
              builder.numInternalNodes );

  aabb_t initBox;
  initBox.lower.x = builder.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.BBox.lower.x + Displace;
  initBox.upper.y = builder.BBox.lower.y + Displace;

  device_queue.submit([&](cl::sycl::handler &cgh) {

    // cgh.parallel_for( cl::sycl::range<1>(nCols*nRows), 
    //                   _query(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints));
    cgh.parallel_for( cl::sycl::range<2>(nCols,nRows), 
                      _query(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints) );

  }).wait_and_throw();

  return;
}


void traverseIterativeCPU(const LBVHoct& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
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

    return;
}


class _queryCPU{
  
  public:
    _queryCPU(const LBVHoct _lbvh, const aabb_t _initBox, const double _Displace, const double _Overlap,
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

            traverseIterativeCPU(lbvh, oCell, pointsCount, tmpId);

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

            traverseIterativeCPU(lbvh, cellBox, pointsCount, idmin);

          }

          count[jj*nCols + ii] = ( minNumPoints <= pointsCount )? idmin : 0u;

        }

      }


      return;
    }
  
  private:
    const LBVHoct  lbvh;
    const aabb_t  initBox;
    const double Overlap;
    const double Displace;
    uint32_t* count;
    const uint32_t nCols;
    // const uint32_t nRows;
    const uint32_t minNumPoints;

};

void stage1query2DCPU(const LBVHoct& lbvh, const whole_t& whole, uint32_t* count, const uint32_t Wsize, 
  const double Overlap, const uint32_t wCols, const uint32_t nCols, const uint32_t nRows, const uint32_t minNumPoints){

  double Displace = round2d(Wsize*(1-Overlap));

  aabb_t initBox;
  initBox.lower.x = whole.lower.x - Wsize + Displace;
  initBox.lower.y = whole.lower.y - Wsize + Displace;
  initBox.upper.x = whole.lower.x + Displace;
  initBox.upper.y = whole.lower.y + Displace;

  // size_t dimChunk = 2;

  tbb::parallel_for( tbb::blocked_range2d<int,int>{ 0, static_cast<int>(nRows),
                                                    static_cast<int>(wCols) ,static_cast<int>(nCols)},
                    _queryCPU(lbvh, initBox, Displace, Overlap, count, nCols, minNumPoints) );

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

  LBVHoct lbvh(
    builder.m_octree,
    builder.m_aabb,
    builder.bintree.ord_point_cloud
  );

  // const size_t dimChunk = 8;

  tbb::parallel_for( tbb::blocked_range2d<int,int>{ 0, static_cast<int>(nRows),
                                                    static_cast<int>(wCols) ,static_cast<int>(nCols)},
                    _queryCPU(lbvh, initBox, Displace, Overlap, count, nCols, minNumPoints) );

  return;
}

void traverseIterative(const LBVHoct& lbvh, aabb_t& queryAABB, uint32_t& numPts, uint32_t& idmin)
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
              real_t min = std::get<0>(pack);
              numPts += std::get<2>(pack);

              if(min < zmin){
                zmin = min;
                idmin = std::get<1>(pack);
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
        // #pragma unroll 8
        for(int i=0; i<8; i++){

          uint32_t child = children[i];

          if ( (overlap & (1u << i)) && lbvh.isLeaf(child) ){

            // const uint32_t l = lbvh.getNumPoints(child);
            // const point_t* p = lbvh.getPointList(child);
            // const uint32_t pointID = lbvh.getObjID(child);

            const auto pack = lbvh.getPointPack(child);
            const point_t* p = std::get<0>(pack);
            // const uint32_t pointID = std::get<1>(pack);
            const uint32_t l = std::get<2>(pack);

            for(int i=0; i<l; i++){
              // point_t p = lbvh.getPoint(child, i);
              if(insideBox(queryAABB, p[i])){
                numPts++;
                if(p[i].z < zmin){
                  zmin = p[i].z;
                  idmin = std::get<1>(pack) + i;
                }
              }
            }

          }

        }

        // __syncthreads();


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

class _queryOct{
  
  public:
    _queryOct(const LBVHoct _lbvh, const aabb_t _initBox, const double _Displace, uint32_t* _count, 
          const uint32_t _nCols, const uint32_t _nRows, const uint32_t _minNumPoints) : 
          lbvh(_lbvh), initBox(_initBox), Displace(_Displace), count(_count),
          nCols(_nCols), nRows(_nRows), minNumPoints(_minNumPoints) {}
    
    void operator()(cl::sycl::id<2> index) const
    {
      int idx = static_cast<int>(index[0]);
      int jdx = static_cast<int>(index[1]);

      if (idx >= nCols || jdx >= nRows)
          return;
      
      aabb_t cellBox;
      cellBox.upper.x = idx*Displace + initBox.upper.x;
      cellBox.lower.x = idx*Displace + initBox.lower.x;
      cellBox.upper.y = jdx*Displace + initBox.upper.y;
      cellBox.lower.y = jdx*Displace + initBox.lower.y;

      uint32_t pointsCount = 0;
      uint32_t idmin = 0xFFFFFFFFu;

      traverseIterative(lbvh, cellBox, pointsCount, idmin);

      count[jdx*nCols + idx] = ( minNumPoints <= pointsCount )? idmin : 0u;

      return;
    }
  
  private:
    const LBVHoct  lbvh;
    const aabb_t  initBox;
    // const double Overlap;
    const double Displace;
    uint32_t* count;
    const uint32_t nCols;
    const uint32_t nRows;
    const uint32_t minNumPoints;

};


sycl::event stage1query2D(Octree_builder& builder, uint32_t* count, uint32_t Wsize, 
  double Overlap, uint32_t nCols, uint32_t nRows, uint32_t minNumPoints,
  cl::sycl::queue device_queue){

  // size_t Ncells = nRows*nCols;

  double Displace = round2d(Wsize*(1-Overlap));

  LBVHoct lbvh(
    builder.m_octree,
    builder.m_aabb,
    builder.bintree.ord_point_cloud
  );

  aabb_t initBox;
  initBox.lower.x = builder.bintree.BBox.lower.x - Wsize + Displace;
  initBox.lower.y = builder.bintree.BBox.lower.y - Wsize + Displace;
  initBox.upper.x = builder.bintree.BBox.lower.x + Displace;
  initBox.upper.y = builder.bintree.BBox.lower.y + Displace;

  return device_queue.submit([&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<2>(nCols,nRows), 
                      _queryOct(lbvh, initBox, Displace, count, nCols, nRows, minNumPoints) );

  });

  // return e;
}


#endif // TRAVERSE_H