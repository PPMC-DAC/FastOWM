
#ifdef CHECK
#include <random>
#endif

namespace sycl_builder {

// __global__ void init_struct( bintree_node* n, aabb_t* bb, const uint32_t num_objects, 
//     const bintree_node default_node, const aabb_t default_aabb)
// {
//   uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
//   uint32_t stride = blockDim.x * gridDim.x;

//   for(uint32_t i = index; i<num_objects; i += stride)
//   {
//     n[i] = default_node;
//     bb[i] = default_aabb;
//   }
//   return;
// }

class _init_struct{
  
  public:
    _init_struct( bintree_node* n, aabb_t* bb, uint32_t* idxs, uint32_t* m, uint32_t* f,
                  const bintree_node _default_node, const aabb_t _default_aabb, 
                  const uint32_t _n, const uint32_t _nin) : 
                  node_list(n), aabb_list(bb), indices(idxs), flags(f), morton(m),
                  default_node(_default_node), default_aabb(_default_aabb), num_nodes(_n),
                  num_internal_nodes(_nin) {}

    void operator()(sycl::id<1> idx) const
    {
      indices[idx] = idx;
      morton[idx] = 0u;

      if(idx >= num_nodes)
        return;

      aabb_list[idx] = default_aabb;
      node_list[idx] = default_node;

      if(idx >= num_internal_nodes) 
        return;
      
      flags[idx] = 0u;
      // // atomic_acc(flags[idx]).store(0u);
      // // sycl::atomic_store(flags[idx], 0u);

      return;
    }
      
  private:
    bintree_node*       node_list;
    aabb_t*             aabb_list;
    uint32_t*           indices;
    uint32_t*           morton;
    uint32_t*           flags;
    const bintree_node  default_node;
    const aabb_t        default_aabb;
    const uint32_t      num_nodes;
    const uint32_t      num_internal_nodes;

};

// __global__ 
// void get_morton_code( const point_t* points, const whole_t BBox, const point_t diffBox, 
//                         uint32_t* morton, uint32_t* indices, const uint32_t num_objects)
// {
//   int index = threadIdx.x + blockIdx.x * blockDim.x;
//   int stride = blockDim.x * gridDim.x;

//   for(int i = index; i<num_objects; i += stride){

//     indices[i] = i;

//     point_t aux = points[i];

//     aux.x -= BBox.lower.x;
//     aux.y -= BBox.lower.y;
//     aux.x /= diffBox.x;
//     aux.y /= diffBox.y;

//     morton[i] = morton2D(aux);

//   }

//   return;
// }

class _get_morton{
  
  public:
    _get_morton( const point_t* p, const whole_t BBox, const point_t diffBox, 
        uint32_t* m) : 
        point_cloud(p), whole(BBox), diff(diffBox), morton(m) {}

    void operator()(sycl::id<1> idx) const
    {
      point_t aux = point_cloud[idx];

      aux.x -= whole.lower.x;
      aux.y -= whole.lower.y;
      aux.x /= diff.x;
      aux.y /= diff.y;

      morton[idx] = morton2D(aux);

      return;
    }
    
    // void operator()(sycl::nd_item<1> it) const
    // {
    //   const uint32_t index = it.get_local_id(0) + it.get_group(0) * it.get_local_range(0);
    //   const uint32_t stride = it.get_local_range(0) * it.get_group_range(0);

    //   for(int i=index; i<numObjects; i+=stride)
    //   {
    //     point_t aux = point_cloud[i];

    //     aux.x -= whole.lower.x;
    //     aux.y -= whole.lower.y;
    //     aux.x /= diff.x;
    //     aux.y /= diff.y;

    //     morton[i] = morton2D(aux);
    //   }

    //   return;
    // }
  
  private:
    const point_t*  point_cloud;
    const whole_t   whole;
    const point_t   diff;
    uint32_t*       morton;

};

// __global__
// void order_points(const point_t* p, point_t* op, const uint32_t* idx, const uint32_t num_objects)
// {
//   const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
//   const uint32_t stride = blockDim.x * gridDim.x;

//   for(uint32_t i = index; i<num_objects; i += stride)
//   {
//     op[i] = p[idx[i]];
//   }

//   return;
// }

class _order_points{
  
  public:
    _order_points( const point_t* p, point_t* op, const uint32_t* idxs) : 
        point_cloud(p), ordered_point_cloud(op), indices(idxs) {}

    void operator()(sycl::id<1> idx) const
    {
      ordered_point_cloud[idx] = point_cloud[indices[idx]];

      return;
    }
    
    // void operator()(sycl::nd_item<1> it) const
    // {
    //   const uint32_t index = it.get_local_id(0) + it.get_group(0) * it.get_local_range(0);
    //   const uint32_t stride = it.get_local_range(0) * it.get_group_range(0);

    //   for(int i=index; i<numObjects; i+=stride)
    //   {
    //     ordered_point_cloud[i] = point_cloud[indices[i]];
    //   }
    //   return;
    // }
  
  private:
    const point_t*  point_cloud;
    point_t*  ordered_point_cloud;
    const uint32_t* indices;

};

// /*con esta puedo evitar usar check_morton, porque a la misma vez inicializo los nodos hoja
// voy calculando los aabb y el morton64 asociado*/
// __global__
// void init_leafs_size(bintree_node* n, aabb_t* bb, const point_t* p,
//   const uint32_t num_objects, const uint32_t num_leafs, const uint32_t leafSize, const aabb_t default_aabb)
// {
//   const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
//   const uint32_t stride = blockDim.x * gridDim.x;

//   for(uint32_t i = index; i<num_leafs; i += stride)
//   {
  
//     const uint32_t start = i*leafSize;
//     const uint32_t end = i*leafSize+leafSize;
//     // const uint32_t end = (start+leafSize <= num_objects)? start+leafSize : num_objects;
//     aabb_t auxbb = default_aabb;

//     real_t min = inf;
//     uint32_t idmin = 0xFFFFFFFFu;

//     /* modifico el índice del nodo hoja par que apunte al
//     índice del primer punto en el vector ordenado de índices,
//     es decir, este NO es el índice del primer punto sino el índice
//     del primer índice en el vector de índices */
//     n[i].object_idx = start;

//     int j=start;
//     for(; j<end && j<num_objects; j++) {
//     // for(; j<end; j++) {
//       // copy the point
//       point_t auxp = p[j];
//       // obtain the new aabb
//       auxbb = merge(auxbb, {{auxp.x,auxp.y},{auxp.x,auxp.y}});
//       // compare min
//       if(auxp.z-min < TOL){
//         min = auxp.z;
//         idmin = j;
//       }
//     }

//     __syncthreads();

//     n[i].numPts = j-start;
//     n[i].min_idx = idmin;
//     bb[i] = auxbb;

//     // __syncthreads();
//   }

//   return;
// }


class _init_leafs{
  
  public:
    _init_leafs( bintree_node* n, aabb_t* bb, point_t* op, const uint32_t num_objects,
        const uint32_t _leafSize, const aabb_t dbb) : 
        node_list(n), aabb_list(bb), ordered_point_cloud(op), numObjects(num_objects),
        leafSize(_leafSize), default_aabb(dbb) {}

    void operator()(sycl::id<1> idx) const
    {
      const uint32_t start = idx*leafSize;
      const uint32_t end = idx*leafSize+leafSize;
      aabb_t auxbb = default_aabb;

      real_t min = inf;
      uint32_t idmin = 0xFFFFFFFFu;

      /* modifico el índice del nodo hoja par que apunte al
      índice del primer punto en el vector ordenado de índices,
      es decir, este NO es el índice del primer punto sino el índice
      del primer índice en el vector de índices */
      node_list[idx].object_idx = start;

      int j=start;
      for(; j<end && j<numObjects; j++) {
        // copy the point
        point_t auxp = ordered_point_cloud[j];
        // obtain the new aabb
        auxbb = merge(auxbb, {{auxp.x,auxp.y},{auxp.x,auxp.y}});
        // compare min
        if(auxp.z-min < TOL){
          min = auxp.z;
          idmin = j;
        }
      }

      node_list[idx].numPts = j-start;
      node_list[idx].min_idx = idmin;
      aabb_list[idx] = auxbb;

      return;
    }
    
//     void operator()(sycl::nd_item<1> it) const
//     {
//       const uint32_t index = it.get_local_id(0) + it.get_group(0) * it.get_local_range(0);
//       const uint32_t stride = it.get_local_range(0) * it.get_group_range(0);

//       for(int i=index; i<numLeafs; i+=stride)
//       {
//         const uint32_t start = i*leafSize;
//         const uint32_t end = i*leafSize+leafSize;
//         aabb_t auxbb = default_aabb;

//         real_t min = inf;
//         uint32_t idmin = 0xFFFFFFFFu;

//         /* modifico el índice del nodo hoja par que apunte al
//         índice del primer punto en el vector ordenado de índices,
//         es decir, este NO es el índice del primer punto sino el índice
//         del primer índice en el vector de índices */
//         node_list[i].object_idx = start;

//         int j=start;
//         for(; j<end && j<numObjects; j++) {
//           // copy the point
//           point_t auxp = ordered_point_cloud[j];
//           // obtain the new aabb
//           auxbb = merge(auxbb, {{auxp.x,auxp.y},{auxp.x,auxp.y}});
//           // compare min
//           if(auxp.z-min < TOL){
//             min = auxp.z;
//             idmin = j;
//           }
//         }

//         node_list[i].numPts = j-start;
//         node_list[i].min_idx = idmin;
//         aabb_list[i] = auxbb;

//       }

//       return;
//     }
  
  private:
    bintree_node*   node_list;
    aabb_t*         aabb_list; 
    const point_t*  ordered_point_cloud;
    const uint32_t  numObjects;
    const uint32_t  leafSize;
    const aabb_t    default_aabb;

};

// __global__
// void init_nodes5( bintree_node* n, const uint32_t num_leafs )
// {
//   const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
//   const uint32_t stride = blockDim.x * gridDim.x;

//   for(int i=index; i<num_leafs-1; i+=stride)
//   {
    
//     uint32_t idx = i;

//     const uint64_t self_code = idx;
//     const int L_delta = (i==0)? FW_S32_MIN : ::__clzll(self_code ^ (idx-1));
//     const int R_delta = ::__clzll(self_code ^ (idx+1));
//     const int d = (R_delta > L_delta) ? 1 : -1;

//     // Compute upper bound for the length of the range

//     const int delta_min = ::min(L_delta, R_delta);
//     uint32_t l_max = 64;
//     int i_tmp = idx + l_max * d;

//     do{

//       l_max<<=1;
//       i_tmp = idx + l_max * d;

//     } while( 0 <= i_tmp && i_tmp < num_leafs && delta_min < ::__clzll(self_code ^ i_tmp) );


//     // Find the other end by binary search
//     uint32_t l = 0;
//     for(uint32_t t= l_max >> 1; t>0; t>>=1)
//     {
//       i_tmp = idx + (l + t) * d;
//       if( 0 <= i_tmp && i_tmp < num_leafs && delta_min < ::__clzll(self_code ^ i_tmp))
//         l += t;
//     }

//     uint32_t jdx = idx + l * d;
//     if(d < 0)
//     {
//         swap(idx, jdx); // make it sure that idx < jdx
//     }

//     const uint64_t first_code = idx;
//     const uint32_t prefix_node = ::__clzll(first_code ^ jdx);

//     // binary search...
//     uint32_t gamma  = idx;
//     uint32_t stride = l;

//     do
//     {
//         stride = (stride + 1) >> 1;
//         const uint32_t middle = gamma + stride;
//         if( middle < jdx && prefix_node < ::__clzll(first_code ^ middle))
//           gamma = middle;

//     } while(stride > 1);

//     // k = gamma

//     uint32_t lidx = (idx == gamma) ? (gamma + num_leafs - 1) : gamma;
//     uint32_t ridx = (jdx == gamma + 1) ? (gamma + num_leafs) : (gamma + 1);

//     n[i].left_idx = lidx;
//     n[i].right_idx = ridx;
//     n[lidx].parent_idx  = i;
//     n[ridx].parent_idx = i;

//     __syncthreads();
//   }

//   return;

// }

class _init_nodes{
  
  public:
    _init_nodes( bintree_node* n, const uint32_t num_leafs) : node_list(n), numLeafs(num_leafs) {}

    void operator()(sycl::id<1> i) const
    {
      uint32_t idx = i;

      const uint64_t self_code = idx;
      const int L_delta = (i==0)? 0 : sycl::clz(self_code ^ (idx-1));
      const int R_delta = sycl::clz(self_code ^ (idx+1));
      const int d = (R_delta > L_delta) ? 1 : -1;

      // Compute upper bound for the length of the range

      const int delta_min = sycl::min(L_delta, R_delta);
      uint32_t l_max = 64;
      int i_tmp = idx + l_max * d;

      do{

        l_max<<=1;
        i_tmp = idx + l_max * d;

      } while( 0 <= i_tmp && i_tmp < numLeafs && delta_min < sycl::clz(self_code ^ i_tmp) );


      // Find the other end by binary search
      uint32_t l = 0;
      for(uint32_t t= l_max >> 1; t>0; t>>=1)
      {
        i_tmp = idx + (l + t) * d;
        if( 0 <= i_tmp && i_tmp < numLeafs && delta_min < sycl::clz(self_code ^ i_tmp))
          l += t;
      }

      uint32_t jdx = idx + l * d;
      if(d < 0)
      {
          swap(idx, jdx); // make it sure that idx < jdx
      }

      const uint64_t first_code = idx;
      const uint32_t prefix_node = sycl::clz(first_code ^ jdx);

      // binary search...
      uint32_t gamma  = idx;
      uint32_t stride = l;

      do
      {
          stride = (stride + 1) >> 1;
          const uint32_t middle = gamma + stride;
          if( middle < jdx && prefix_node < sycl::clz(first_code ^ middle))
            gamma = middle;

      } while(stride > 1);

      // k = gamma

      uint32_t lidx = (idx == gamma) ? (gamma + numLeafs - 1) : gamma;
      uint32_t ridx = (jdx == gamma + 1) ? (gamma + numLeafs) : (gamma + 1);

      node_list[i].left_idx = lidx;
      node_list[i].right_idx = ridx;
      node_list[lidx].parent_idx  = i;
      node_list[ridx].parent_idx = i;

      return;
    }
    
//     void operator()(sycl::nd_item<1> it) const
//     {
//       const uint32_t index = it.get_local_id(0) + it.get_group(0) * it.get_local_range(0);
//       const uint32_t stride = it.get_local_range(0) * it.get_group_range(0);
//       const uint32_t num_internal_nodes = numLeafs-1;

//       for(int i=index; i<num_internal_nodes; i+=stride)
//       {
        
//         uint32_t idx = i;

//         const uint64_t self_code = idx;
//         const int L_delta = (i==0)? 0 : sycl::clz(self_code ^ (idx-1));
//         const int R_delta = sycl::clz(self_code ^ (idx+1));
//         const int d = (R_delta > L_delta) ? 1 : -1;

//         // Compute upper bound for the length of the range

//         const int delta_min = sycl::min(L_delta, R_delta);
//         uint32_t l_max = 64;
//         int i_tmp = idx + l_max * d;

//         do{

//           l_max<<=1;
//           i_tmp = idx + l_max * d;

//         } while( 0 <= i_tmp && i_tmp < numLeafs && delta_min < sycl::clz(self_code ^ i_tmp) );


//         // Find the other end by binary search
//         uint32_t l = 0;
//         for(uint32_t t= l_max >> 1; t>0; t>>=1)
//         {
//           i_tmp = idx + (l + t) * d;
//           if( 0 <= i_tmp && i_tmp < numLeafs && delta_min < sycl::clz(self_code ^ i_tmp))
//             l += t;
//         }

//         uint32_t jdx = idx + l * d;
//         if(d < 0)
//         {
//             swap(idx, jdx); // make it sure that idx < jdx
//         }

//         const uint64_t first_code = idx;
//         const uint32_t prefix_node = sycl::clz(first_code ^ jdx);

//         // binary search...
//         uint32_t gamma  = idx;
//         uint32_t stride = l;

//         do
//         {
//             stride = (stride + 1) >> 1;
//             const uint32_t middle = gamma + stride;
//             if( middle < jdx && prefix_node < sycl::clz(first_code ^ middle))
//               gamma = middle;

//         } while(stride > 1);

//         // k = gamma

//         uint32_t lidx = (idx == gamma) ? (gamma + numLeafs - 1) : gamma;
//         uint32_t ridx = (jdx == gamma + 1) ? (gamma + numLeafs) : (gamma + 1);

//         node_list[i].left_idx = lidx;
//         node_list[i].right_idx = ridx;
//         node_list[lidx].parent_idx  = i;
//         node_list[ridx].parent_idx = i;

//       }

//       return;
//     }
  
  private:
    bintree_node*  node_list;
    const uint32_t numLeafs;

};


// /*
// Igual que el caso del kernel anterior, en este tampoco puedo lanzar un thread por
// cada nodo hoja y tengo que utilizar un stride para que cada thread haga más trabajo
// */
// __global__
// void create_aabb_size( bintree_node* n, aabb_t* bb, uint32_t* flags, const point_t* p, 
//   const uint32_t num_internal_nodes, const uint32_t num_nodes)
// {
//   uint32_t index = num_internal_nodes + threadIdx.x + blockIdx.x * blockDim.x;
//   uint32_t stride = blockDim.x * gridDim.x;

//   for(uint32_t i = index; i < num_nodes; i += stride)
//   {
//     uint32_t parent = n[i].parent_idx;

//     while(parent != 0xFFFFFFFFu && atomicCAS(flags + parent, 0, 1)) 
//     {
//       // here, the flag has already been 1. it means that this
//       // thread is the 2nd thread. merge AABB of both childlen.

//       const uint32_t lidx = n[parent].left_idx;
//       const uint32_t ridx = n[parent].right_idx;
//       const aabb_t lbox = bb[lidx];
//       const aabb_t rbox = bb[ridx];
//       const uint32_t lmin = n[lidx].min_idx;
//       const uint32_t rmin = n[ridx].min_idx;
      
//       // compare mins
//       if(p[lmin].z - p[rmin].z < TOL)
//         n[parent].min_idx = lmin;
//       else
//         n[parent].min_idx = rmin;

//       // count points
//       n[parent].numPts = n[lidx].numPts + n[ridx].numPts;
//       // merge aabbs
//       bb[parent] = merge(lbox, rbox);

//       // look the next parent...
//       parent = n[parent].parent_idx;
//     }

//     __syncthreads();
//   }

//   return;
// }

class _create_aabbs{
  
  public:
    _create_aabbs( bintree_node* n, aabb_t* bb, uint32_t* f, const point_t* op,
                    uint32_t nn, uint32_t nin) : 
        node_list(n), aabb_list(bb), flags(f), ordered_point_cloud(op),
        numNodes(nn), numInternalNodes(nin) {}

    // void operator()(sycl::id<1> idx) const
    // {

    //   uint32_t parent = node_list[idx].parent_idx;

    //   while(parent != 0xFFFFFFFFu && atomic_acc(flags[parent]).exchange(1u)) 
    //   {

    //     const uint32_t lidx = node_list[parent].left_idx;
    //     const uint32_t ridx = node_list[parent].right_idx;
    //     const aabb_t lbox = aabb_list[lidx];
    //     const aabb_t rbox = aabb_list[ridx];
    //     const uint32_t lmin = node_list[lidx].min_idx;
    //     const uint32_t rmin = node_list[ridx].min_idx;

    //     // compare mins
    //     if(ordered_point_cloud[lmin].z - ordered_point_cloud[rmin].z < TOL)
    //       node_list[parent].min_idx = lmin;
    //     else
    //       node_list[parent].min_idx = rmin;

    //     // count points
    //     node_list[parent].numPts = node_list[lidx].numPts + node_list[ridx].numPts;
    //     // merge aabbs
    //     aabb_list[parent] = merge(lbox, rbox);

    //     // look the next parent...
    //     parent = node_list[parent].parent_idx;
    //   }

    //   return;

    // }// operator() id
    
    void operator()(sycl::nd_item<1> it) const
    {

      uint32_t index = numInternalNodes + it.get_local_id(0) + it.get_group(0) * it.get_local_range(0);
      uint32_t stride = it.get_local_range(0) * it.get_group_range(0);

      for(uint32_t i = index; i < numNodes; i += stride)
      {
        uint32_t parent = node_list[i].parent_idx;

        while(parent != 0xFFFFFFFFu && atomic_acc(flags[parent]).exchange(1u)) 
        {
          // here, the flag has already been 1. it means that this
          // thread is the 2nd thread. merge AABB of both childlen.

          const uint32_t lidx = node_list[parent].left_idx;
          const uint32_t ridx = node_list[parent].right_idx;
          const aabb_t lbox = aabb_list[lidx];
          const aabb_t rbox = aabb_list[ridx];
          const uint32_t lmin = node_list[lidx].min_idx;
          const uint32_t rmin = node_list[ridx].min_idx;
          
          // compare mins
          if(ordered_point_cloud[lmin].z - ordered_point_cloud[rmin].z < TOL)
            node_list[parent].min_idx = lmin;
          else
            node_list[parent].min_idx = rmin;

          // count points
          node_list[parent].numPts = node_list[lidx].numPts + node_list[ridx].numPts;
          // merge aabbs
          aabb_list[parent] = merge(lbox, rbox);

          // look the next parent...
          parent = node_list[parent].parent_idx;
        }

        // it.barrier();

      }

      return;
      
    } // operator() nd_item
  
  private:
    bintree_node* node_list;
    aabb_t* aabb_list;
    uint32_t* flags;
    const point_t*  ordered_point_cloud;
    uint32_t numNodes;
    uint32_t numInternalNodes;

};


void* mallocWrap(const size_t& size, sycl::queue device_queue)
{
#ifdef SHARED
  void *ptr = malloc_shared(size, device_queue);
#elif DEVICE
  void *ptr = malloc_device(size, device_queue);
#else
  void *ptr = malloc_host(size, device_queue);
#endif

  if (ptr)
      return ptr;
  else
      throw std::bad_alloc{};
}


} // namespace sycl_builder


SYCL_builder::SYCL_builder( std::string inputTXT, const uint32_t chunk, sycl::queue& q) : 
  leafSize(chunk), device_queue(q)
{
#ifndef CHECK
  readHeader(inputTXT, BBox, numObjects);
#else
  numObjects = 32;

  std::mt19937 mt(123456789);
  std::uniform_real_distribution<real_t> uni(0.0, 10.0);

  BBox.lower.x = 0.0;
  BBox.lower.y = 0.0;
  BBox.upper.x = 10.0;
  BBox.upper.y = 10.0;
#endif

  numLeafs = (uint32_t)((numObjects-1)/leafSize) + 1u;
  numInternalNodes = numLeafs - 1u;
  numNodes = 2*numLeafs - 1u;
  
  diffBox.x = BBox.upper.x - BBox.lower.x;
  diffBox.y = BBox.upper.y - BBox.lower.y;

  num_groups = 8;
  items_per_group = 64;
  grid_dim = num_groups * items_per_group;

  point_cloud = (point_t*)sycl_builder::mallocWrap(numObjects*sizeof(point_t), device_queue);
  ord_point_cloud = (point_t*)sycl_builder::mallocWrap(numObjects*sizeof(point_t), device_queue);

  morton = (uint32_t*)sycl_builder::mallocWrap(numObjects*sizeof(uint32_t), device_queue);
  indices = (uint32_t*)sycl_builder::mallocWrap(numObjects*sizeof(uint32_t), device_queue);

  node_list = (bintree_node*)sycl_builder::mallocWrap(numNodes * sizeof(bintree_node), device_queue);
  aabb_list = (aabb_t*)sycl_builder::mallocWrap(numNodes * sizeof(aabb_t), device_queue);

  flags = (uint32_t*)sycl_builder::mallocWrap(numInternalNodes*sizeof(uint32_t), device_queue);
  // flags = (uint32_t*)sycl::malloc_device(numInternalNodes*sizeof(uint32_t), device_queue);

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, NUM_PROCS);
  
#ifdef DEBUG
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif

#ifdef DEVICE
  point_t* point_cloud_h = (point_t*)malloc(numObjects*sizeof(point_t));
#endif

#ifndef CHECK

#ifdef DEVICE
  map_file_atof_tbb_th(inputTXT, point_cloud_h, numObjects);
#else
  map_file_atof_tbb_th(inputTXT, point_cloud, numObjects);
#endif

#else

#ifdef DEVICE
  for(int i=0; i<numObjects; i++){
    point_cloud_h[i].x = (int)uni(mt);
    point_cloud_h[i].y = (int)uni(mt);
    point_cloud_h[i].z = (int)uni(mt);
  }
  /*este tiene que ser siempre el mínimo en root si todo ha ido bien*/
  point_cloud_h[9].z = -1.0;
#else
  for(int i=0; i<numObjects; i++){
    point_cloud[i].x = (int)uni(mt);
    point_cloud[i].y = (int)uni(mt);
    point_cloud[i].z = (int)uni(mt);
  }
  /*este tiene que ser siempre el mínimo en root si todo ha ido bien*/
  point_cloud[9].z = -1.0;
#endif

#endif

#ifdef DEVICE
  device_queue.memcpy(point_cloud, point_cloud_h, numObjects*sizeof(point_t));
  device_queue.wait_and_throw();

  free(point_cloud_h);
#endif

#ifdef DEBUG
  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << "Load time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  std::cout << "INIT "; 
  i_start = tempo_t::now();
#endif

  reset();

#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << "time elapsed: " << mytime << " ms\n";
#endif

} // constructor


SYCL_builder::~SYCL_builder()
{

  free(point_cloud, device_queue);
  free(ord_point_cloud, device_queue);
  free(morton, device_queue);
  free(indices, device_queue);
  free(node_list, device_queue);
  free(aabb_list, device_queue);
  free(flags, device_queue);

} // destructor



void SYCL_builder::build()
{

#ifdef DEBUG
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif


  device_queue.submit([&](sycl::handler &cgh) {

    cgh.parallel_for(sycl::range<1>(numObjects), 
                      sycl_builder::_get_morton(point_cloud, BBox, diffBox, morton));

  }).wait_and_throw();

	
#ifdef DEBUG
  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " MORTON time elapsed: " << mytime << " ms\n";

  i_start = tempo_t::now();
#endif

#ifdef CHECK
  for(int i = 0; i<numObjects; i++){
    std::cout << morton[i] << " ";
  }
  std::cout << "\n";
  for(int i = 0; i<numObjects; i++){
    std::cout << indices[i] << " ";
  }
  std::cout << "\n";
#endif

  auto itt = oneapi::dpl::make_zip_iterator(morton, indices);

#if CPU | FEMU
  oneapi::dpl::stable_sort(	oneapi::dpl::execution::par_unseq,
                      itt,
                      itt + numObjects,
                      [](auto a, auto b){ 
                        return std::get<0>(a) < std::get<0>(b);
                      });
#elif GPU

  auto policy = oneapi::dpl::execution::make_device_policy(device_queue);

  oneapi::dpl::stable_sort(	policy,
                      itt,
                      itt + numObjects,
                      [](auto a, auto b){ 
                        return std::get<0>(a) < std::get<0>(b);
                      });
#endif

#ifdef CHECK
  for(int i = 0; i<numObjects; i++){
    std::cout << morton[i] << " ";
  }
  std::cout << "\n";
  for(int i = 0; i<numObjects; i++){
    std::cout << indices[i] << " ";
  }
  std::cout << "\n";
#endif
	
#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " SORT time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif


  sycl::event e1 = device_queue.submit([&](sycl::handler &cgh) {

    cgh.parallel_for(sycl::range<1>(numObjects), 
                      sycl_builder::_order_points(point_cloud, ord_point_cloud, indices));

  });


#ifdef DEBUG

  e1.wait();

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " ORDER time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif

  sycl::event e2 = device_queue.submit([&](sycl::handler &cgh) {

    cgh.depends_on(e1);

    cgh.parallel_for(sycl::range<1>(numLeafs), 
                      sycl_builder::_init_leafs(&node_list[numInternalNodes], &aabb_list[numInternalNodes], 
                                  ord_point_cloud, numObjects, leafSize, default_aabb));

  });

#ifdef DEBUG

  e2.wait();

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " LEAFS and AABBs time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif
																		  

  sycl::event e3 = device_queue.submit([&](sycl::handler &cgh) {
    
    cgh.parallel_for(sycl::range<1>(numInternalNodes), 
                      sycl_builder::_init_nodes(node_list, numLeafs));

  });


#ifdef DEBUG

  e3.wait();

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " INTERNAL NODES time elapsed: " << mytime << " ms\n";
#endif

#ifdef CHECK
  for(int i = 0; i<numNodes; i++){
    std::cout << node_list[i].parent_idx << ", " << node_list[i].left_idx << ", ";
    std::cout << node_list[i].right_idx << ", " << node_list[i].object_idx << ", ";
    std::cout << node_list[i].numPts << ", " << node_list[i].min_idx << "\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<numNodes; i++){
      std::cout << aabb_list[i].upper.x << "," << aabb_list[i].upper.y << " ";
      std::cout << aabb_list[i].lower.x << "," << aabb_list[i].lower.y << "\n";
    }
  std::cout << std::endl;
  for(int i = 0; i<numInternalNodes; i++){
    std::cout << flags[i] << " ";
  }
  std::cout << "\n\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif
									
  device_queue.submit([&](sycl::handler &cgh) {

    cgh.depends_on({e2,e3});

    cgh.parallel_for(sycl::nd_range<1>(grid_dim, items_per_group), 
                    sycl_builder::_create_aabbs(node_list, aabb_list, flags, ord_point_cloud, numNodes, numInternalNodes));

  });

  device_queue.wait_and_throw();


#ifdef CHECK
  for(int i = 0; i<numNodes; i++){
    std::cout << node_list[i].parent_idx << ", " << node_list[i].left_idx << ", ";
    std::cout << node_list[i].right_idx << ", " << node_list[i].object_idx << ", ";
    std::cout << node_list[i].numPts << ", " << node_list[i].min_idx << "(";
    if(node_list[i].min_idx != 0xFFFFFFFFu)
      std::cout << ord_point_cloud[node_list[i].min_idx].z;
    std::cout << ")\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<numNodes; i++){
    std::cout << aabb_list[i].upper.x << "," << aabb_list[i].upper.y << " ";
    std::cout << aabb_list[i].lower.x << "," << aabb_list[i].lower.y << "\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<numInternalNodes; i++){
    std::cout << flags[i] << " ";
  }
  std::cout << "\n\n";
#endif

	
#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " AABBs time elapsed: " << mytime << " ms\n";
#endif

  return;

} // build()

void SYCL_builder::reset()
{

  device_queue.submit([&](sycl::handler &cgh) {

    cgh.parallel_for( sycl::range<1>(numObjects), 
                      sycl_builder::_init_struct(node_list, aabb_list, indices, morton, flags,
                                      default_node, default_aabb, numNodes, numInternalNodes));

  }).wait_and_throw();
  
} // reset()