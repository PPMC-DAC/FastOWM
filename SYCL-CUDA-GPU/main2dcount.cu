
// #include "include/deftypes.cuh"
// #include "search/minsearch.cuh"
#include "minsearch.hpp"

tree array_nodes = NULL;
tree2 array_nodes2 = NULL;
point_t* array_all_points = NULL;
int copied_nodes = 0;
int num_pts_copiados = 0;

size_t num_objects;
uint32_t r_num_objects;
uint32_t  r_num_internal_nodes;
uint32_t r_num_nodes;



__device__ __host__
inline aabb_t merge(const aabb_t& lhs, const aabb_t& rhs) noexcept
{
    aabb_t merged;
    merged.upper.x = ::fmax(lhs.upper.x, rhs.upper.x);
    merged.upper.y = ::fmax(lhs.upper.y, rhs.upper.y);
    merged.lower.x = ::fmin(lhs.lower.x, rhs.lower.x);
    merged.lower.y = ::fmin(lhs.lower.y, rhs.lower.y);
    return merged;
}

struct aabb_getter
{
  /* aqui destaca que el limite superior e inferior que le da a las cajas es el mismo,
  es decir, la caja es un solo punto */
  __device__
  aabb_t operator()(const real_v f) const noexcept
  {
      aabb_t retval;
      retval.upper = f;
      retval.lower = f;
      return retval;
  }
};


int read_pointsC(std::string file_name, std::vector<point_t>& point_cloud)
{
  FILE* fileLAS;

  if((fileLAS = fopen(file_name.c_str(),"r")) == NULL){
    printf("Unable to open file!\n");
    return -1;
  }

  for(auto& p : point_cloud){
    //Obtengo los datos id X Y Z
    if(fscanf(fileLAS, "%lf %lf %lf", &p.x, &p.y, &p.z) < 3){
      printf("Imposible to obtain values\n");
      return -1;
    }
    while(fgetc(fileLAS)!='\n');
  }

  //Ya no necesito mas el fichero
  if(fclose(fileLAS)){
    printf("Cannot close the file\n");
    return -1;
  }

  return 0;
}

__device__
uint32_t SeparateBy1(uint32_t x) {
    x &= 0x0000ffffu;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x <<  8)) & 0x00ff00ffu; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x <<  4)) & 0x0f0f0f0fu; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x <<  2)) & 0x33333333u; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x <<  1)) & 0x55555555u; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    return x;
}

__device__
inline uint32_t morton2D(point_t xyz , const real_t resolution = 1024.0f ) noexcept
{
  xyz.x = ::fmin(xyz.x * resolution, resolution - 1.0);
  xyz.y = ::fmin(xyz.y * resolution, resolution - 1.0);

  const uint32_t xx = SeparateBy1(static_cast<uint32_t>(xyz.x));
  const uint32_t yy = SeparateBy1(static_cast<uint32_t>(xyz.y));
  return xx * 2 + yy;
}

__device__
inline uint32_t morton2D(real_v& xy , const real_t resolution = 1024.0f ) noexcept
{
  // xy.x = ::fmin(xy.x * resolution, resolution - 1.0);
  // xy.y = ::fmin(xy.y * resolution, resolution - 1.0);

  const uint32_t xx = SeparateBy1(static_cast<uint32_t>(xy.x));
  const uint32_t yy = SeparateBy1(static_cast<uint32_t>(xy.y));
  return xx * 2 + yy;
}


__global__ 
void get_morton_code(const point_t* points, const aabb_t BBox, const real_v diffBox, 
                        uint32_t* morton, uint32_t* indices, const size_t N)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for(int i = index; i<N; i += stride){

    indices[i] = i;

    point_t aux = points[i];

#ifndef DEBUG
    aux.x -= BBox.lower.x;
    aux.y -= BBox.lower.y;
    aux.x /= diffBox.x;
    aux.y /= diffBox.y;
#endif

    morton[i] = morton2D(aux);

  }

  return;
}

__global__
void check_morton( const uint32_t* m, const uint32_t* idx, uint64_t* morton64, const size_t N)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  // uint64_t m64;

  for(int i = index; i<N; i += stride)
  {
    uint64_t m64 = m[i];
    m64 <<= 32;
    m64 |= idx[i];
    morton64[i] = m64;
  }

  return;
}

__global__
void check_morton_size( const uint32_t* m, const uint32_t* idx, uint64_t* morton64, 
  const size_t r_objects, const int chunkDim)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  // uint64_t m64;

  for(int i = index; i<r_objects; i += stride)
  {
    uint64_t m64 = m[i*chunkDim];
    m64 <<= 32;
    m64 |= idx[i*chunkDim];
    morton64[i] = m64;
  }

  return;
}

__global__
void init_leafs(node_t* n, aabb_t* bb, const point_t* p, const uint32_t* idx, const size_t N)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for(int i = index; i<N; i += stride)
  {
    uint32_t id = idx[i]; // obtengo el índice del punto en esa posición

    point_t auxp = p[id]; // obtengo el punto

    n[i].object_idx = id; // modifico el nodo
    // auxn.object_idx = id; // modifico el nodo

    bb[i].upper = {auxp.x,auxp.y}; // modifico aabb
    bb[i].lower = {auxp.x,auxp.y};

    // n[i] = auxn;
  }

  return;
}


__global__
void init_leafs_size(node_t* n, aabb_t* bb, const point_t* p, const uint32_t* idx,
  const int num_objects, const int r_objects, const int chunkDim, const aabb_t default_aabb)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  point_t auxp;

  for(int i = index; i<r_objects; i += stride)
  {
    uint32_t start = i*chunkDim;
    uint32_t end = i*chunkDim+chunkDim;
    aabb_t auxbb = default_aabb;

    real_t min =  NPP_MAXABS_64F;
    uint32_t idmin = 0xFFFFFFFFu;
  
    /* modifico el índice del nodo hoja par que apunte al primer punto */
    n[i].object_idx = start;

    // point_t auxp = p[idx[start]]; // obtengo el punto
    // aabb_t auxbb = {{auxp.x,auxp.y},{auxp.x,auxp.y}};
    // if(auxp.z-min < TOL){
    //   min = auxp.z;
    //   idmin = idx[start];
    // }

    int j=start;
    for(; j<end && j<num_objects; j++) {
      // copy the point
      auxp = p[idx[j]];
      // obtain the new aabb
      auxbb = merge(auxbb, {{auxp.x,auxp.y},{auxp.x,auxp.y}});
      // compare min
      if(auxp.z-min < TOL){
        min = auxp.z;
        idmin = idx[j];
      }
    }

    // if(idmin == 0xFFFFFFFFu)
    //   printf("Min set ERROR");

    n[i].numPts = j-start;
    n[i].min_idx = idmin;
    bb[i] = auxbb;
  }

  return;
}

__global__
void init_leafs_size(node_t* n, aabb_t* bb, const point_t* p, point_t* p_out, const uint32_t* idx,
  const int num_objects, const int r_objects, const int chunkDim, const aabb_t default_aabb)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  point_t auxp;

  for(int i = index; i<r_objects; i += stride)
  {
    uint32_t start = i*chunkDim;
    uint32_t end = i*chunkDim+chunkDim;
    aabb_t auxbb = default_aabb;

    real_t min =  NPP_MAXABS_64F;
    uint32_t idmin = 0xFFFFFFFFu;
  
    /* modifico el índice del nodo hoja par que apunte al primer punto */
    n[i].object_idx = start;

    // point_t auxp = p[idx[start]]; // obtengo el punto
    // aabb_t auxbb = {{auxp.x,auxp.y},{auxp.x,auxp.y}};
    // if(auxp.z-min < TOL){
    //   min = auxp.z;
    //   idmin = idx[start];
    // }

    int j=start;
    for(; j<end && j<num_objects; j++) {
      // copy the point
      auxp = p[idx[j]];
      // obtain the new aabb
      auxbb = merge(auxbb, {{auxp.x,auxp.y},{auxp.x,auxp.y}});

      p_out[j] = auxp;
      // compare min
      if(auxp.z-min < TOL){
        min = auxp.z;
        idmin = idx[j];
      }
    }

    // if(idmin == 0xFFFFFFFFu)
    //   printf("Min set ERROR");

    n[i].numPts = j-start;
    n[i].min_idx = idmin;
    bb[i] = auxbb;
  }

  return;
}

__device__ __host__
inline real_v centroid(const aabb_t& box) noexcept
{
    real_v c;
    c.x = (box.upper.x + box.lower.x) * 0.5;
    c.y = (box.upper.y + box.lower.y) * 0.5;
    // c.z = (box.upper.z + box.lower.z) * 0.5;
    return c;
}

/*con esta puedo evitar usar check_morton, porque a la misma vez inicializo los nodos hoja
voy calculando los aabb y el morton64 asociado*/
__global__
void init_leafs_size2(node_t* n, aabb_t* bb, const point_t* p, const uint32_t* idx, uint64_t* m64,
  const int num_objects, const int num_leafs, const int chunkDim, const aabb_t default_aabb,
  const aabb_t BBox, const real_v diffBox)
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  const uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_leafs; i += stride)
  {
    uint32_t start = i*chunkDim;
    uint32_t end = i*chunkDim+chunkDim;
    aabb_t auxbb = default_aabb;

    real_t min =  NPP_MAXABS_64F;
    uint32_t idmin = 0xFFFFFFFFu;
  
    /* modifico el índice del nodo hoja par que apunte al primer punto */
    n[i].object_idx = start;

    // point_t auxp = p[idx[start]]; // obtengo el punto
    // aabb_t auxbb = {{auxp.x,auxp.y},{auxp.x,auxp.y}};
    // if(auxp.z-min < TOL){
    //   min = auxp.z;
    //   idmin = idx[start];
    // }

    int j=start;
    for(; j<end && j<num_objects; j++) {
      // copy the point
      point_t auxp = p[idx[j]];
      // obtain the new aabb
      auxbb = merge(auxbb, {{auxp.x,auxp.y},{auxp.x,auxp.y}});
      // compare min
      if(auxp.z-min < TOL){
        min = auxp.z;
        idmin = idx[j];
      }
    }

    /*Si no estoy utilizando una resolución en la creación de 
    los códigos Morton, hacer todo esto de a continución es lo mismo
    que simplemente guardar el índice. Con la normalización obtengo
    un punto en (0,0)-(1,1). Al truncar es muy probable que el valor
    de la coordenada sea en la mayoría de los casos 0. Tras esto,
    toda operación morton2D() va a dar como resultado un aux64=0 al cual
    le estamos colocando el índice en la parte final, es decir, solo
    nos quedamos con el mínimo en la inmensa mayoría de los casos, aunque
    puede que haya algun punto que sí tuviera un 1 en los primeros 32 bits,
    pero es despreciable*/

    // real_v c = centroid(auxbb);

    // c.x -= BBox.lower.x;
    // c.y -= BBox.lower.y;
    // c.x /= diffBox.x;
    // c.y /= diffBox.y;

    // uint64_t aux64 = morton2D(c);
    // aux64 <<= 32;
    // aux64 |= i;
    // m64[i] = aux64;

    m64[i] = uint64_t(i);

    // if(idmin == 0xFFFFFFFFu)
    //   printf("Min set ERROR");

    n[i].numPts = j-start;
    n[i].min_idx = idmin;
    bb[i] = auxbb;
  }

  return;
}


// __global__
// void init_leafs2(node_t* n, aabb_t* bb, real_v* p, uint32_t* idx, 
//   const size_t nin, const size_t NN, const aabb_t default_aabb)
// {
//   int index = threadIdx.x + blockIdx.x * blockDim.x;
//   int stride = blockDim.x * gridDim.x;

//   node_t auxn;
//   auxn.parent_idx = 0xFFFFFFFF;
//   auxn.left_idx   = 0xFFFFFFFF;
//   auxn.right_idx  = 0xFFFFFFFF;
//   auxn.object_idx = 0xFFFFFFFF;

//   for(int i = index; i<nin; i += stride)
//   {
//     n[i] = auxn;
//     bb[i] = default_aabb;
//   }

//   for(int i = index + nin; i<NN; i += stride)
//   {
//     uint32_t id = idx[i]; // obtengo el índice del punto en esa posición

//     real_v auxp = p[id]; // obtengo el punto

//     // creo un BB que es solo el punto
//     aabb_t auxbb;
//     auxbb.upper = auxp;
//     auxbb.lower = auxp;

//     auxn.object_idx = id; // seteo el índice del nodo

//     bb[i] = auxbb;
//     n[i] = auxn;
//   }

//   return;
// }

// __global__
// void init_indices(node_t* n, uint32_t* idx, const size_t N)
// {
//   int index = threadIdx.x + blockIdx.x * blockDim.x;
//   int stride = blockDim.x * gridDim.x;

//   for(int i = index; i<N; i += stride)
//   {
//     n[i].object_idx = idx[i];
//   }

//   return;
// }


// /* Este pequeño código sirve para demostrar que NVCC permite compilar
// kernels recursivos */
// // __device__ 
// // void recursivo(int id)
// // {
// //   if(id==5){

// //   }else{
// //     recursivo(id+1);
// //   }

// //   return;
// // }
// // __global__
// // void launch_recursivo()
// // {
// //   recursivo(0);
// //   return;
// // }

__global__
void init_nodes( node_t* n, const uint64_t* m64, const size_t r_num_objects )
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  const uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<r_num_objects-1; i += stride)
  {
    const uint2 ij  = lbvh::detail::determine_range(m64, r_num_objects, i);
    // const uint32_t gamma = lbvh::detail::find_split(m64, r_num_objects, ij.x, ij.y);
    uint32_t lidx = lbvh::detail::find_split(m64, r_num_objects, ij.x, ij.y);

    uint32_t ridx = lidx + 1;

    if(min(ij.x, ij.y) == lidx)
    {
        lidx += r_num_objects - 1;
    }
    if(max(ij.x, ij.y) == ridx)
    {
        ridx += r_num_objects - 1;
    }
    // node_listS[idx] = aux;

    n[i].left_idx = lidx;
    n[i].right_idx = ridx;
    n[lidx].parent_idx  = i;
    n[ridx].parent_idx = i;
    return;
}
  return;

}

__device__
inline void swap(uint32_t& left, uint32_t& right)
{
  uint32_t aux = left;
  left = right;
  right = aux;
}

__global__
void init_nodes2( node_t* n, const uint64_t* node_code, const uint32_t num_leafs )
{
  const uint32_t i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i >= num_leafs-1)
    return;
    
  // const uint2 ij = build_internal(node_code, num_leafs, idx);

  uint32_t idx = i;

  const uint64_t self_code = node_code[idx];
  const int L_delta = (i==0)? FW_S32_MIN : ::__clzll(self_code ^ node_code[idx-1]);
  const int R_delta = ::__clzll(self_code ^ node_code[idx+1]);
  const int d = (R_delta > L_delta) ? 1 : -1;

  // Compute upper bound for the length of the range

  // const int delta_min = ::__clzll(self_code ^ node_code[idx-d]);
  const int delta_min = min(L_delta,R_delta);
  uint32_t l_max = 64;
  int i_tmp = idx + l_max * d;

  do{

    l_max<<=1;
    i_tmp = idx + l_max * d;

  } while( 0 <= i_tmp && i_tmp < num_leafs && delta_min < ::__clzll(self_code ^ node_code[i_tmp]) );


  // Find the other end by binary search
  uint32_t l = 0;
  for(uint32_t t= l_max >> 1; t>0; t>>=1)
  {
    i_tmp = idx + (l + t) * d;
    if( 0 <= i_tmp && i_tmp < num_leafs && delta_min < ::__clzll(self_code ^ node_code[i_tmp]))
      l += t;
  }

  uint32_t jdx = idx + l * d;
  if(d < 0)
  {
      swap(idx, jdx); // make it sure that idx < jdx
  }

  const uint64_t first_code = node_code[idx];
  // const uint64_t last_code  = node_code[jdx];
  const uint32_t prefix_node = ::__clzll(first_code ^ node_code[jdx]);

  // binary search...
  uint32_t gamma  = idx;
  uint32_t stride = l;

  do
  {
      stride = (stride + 1) >> 1;
      const uint32_t middle = gamma + stride;
      if( middle < jdx && prefix_node < ::__clzll(first_code ^ node_code[middle]))
        gamma = middle;

  } while(stride > 1);

  // k = gamma

#ifdef DEBUG
  printf("ID: %u, L_d: %d, R_d: %d, d: %d, range: [%u,%u], gamma: %u\n",
          i, L_delta, R_delta, d, idx, jdx, gamma);
#endif


  uint32_t lidx = (idx == gamma) ? (gamma + num_leafs - 1) : gamma;
  uint32_t ridx = (jdx == gamma + 1) ? (gamma + num_leafs) : (gamma + 1);

  n[i].left_idx = lidx;
  n[i].right_idx = ridx;
  n[lidx].parent_idx  = i;
  n[ridx].parent_idx = i;


  return;

}

__global__
void create_aabb( const node_t* n, aabb_t* bb, int* flags, 
      const size_t num_internal_nodes, const size_t num_nodes)
{
  uint32_t index = num_internal_nodes + threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i < num_nodes; i += stride)
  {
    uint32_t parent = n[i].parent_idx;
    while(parent != 0xFFFFFFFF) // means idx == 0
    {
      const int old = atomicCAS(flags + parent, 0, 1);
      if(old == 0)
      {
          // this is the first thread entered here.
          // wait the other thread from the other child node.
          return;
      }
      assert(old == 1);
      // here, the flag has already been 1. it means that this
      // thread is the 2nd thread. merge AABB of both childlen.

      const uint32_t lidx = n[parent].left_idx;
      const uint32_t ridx = n[parent].right_idx;
      const aabb_t lbox = bb[lidx];
      const aabb_t rbox = bb[ridx];
      bb[parent] = merge(lbox, rbox);

      // look the next parent...
      parent = n[parent].parent_idx;
    }

  }

  return;
}

__global__
void create_aabb_size( node_t* n, aabb_t* bb, int* flags, const point_t* p, 
      const size_t num_internal_nodes, const size_t num_nodes)
{
  uint32_t index = num_internal_nodes + threadIdx.x + blockIdx.x * blockDim.x;
  // uint32_t stride = blockDim.x * gridDim.x;
  // printf("index %u, threadIdx, %d, blockIdx %d\n", index, threadIdx.x, blockIdx.x);

  if(index < num_nodes)
  {
    uint32_t parent = n[index].parent_idx;
    while(parent != 0xFFFFFFFF) // means idx == 0
    {
      const int old = atomicCAS(flags + parent, 0, 1);
      if(old == 0)
      {
          // this is the first thread entered here.
          // wait the other thread from the other child node.
          return;
      }
      assert(old == 1);
      // here, the flag has already been 1. it means that this
      // thread is the 2nd thread. merge AABB of both childlen.

      const uint32_t lidx = n[parent].left_idx;
      const uint32_t ridx = n[parent].right_idx;
      const aabb_t lbox = bb[lidx];
      const aabb_t rbox = bb[ridx];
      const uint32_t lmin = n[lidx].min_idx;
      const uint32_t rmin = n[ridx].min_idx;
      // compare mins
#ifdef DEBUG
      if(lmin == 0xFFFFFFFFu){
        uint32_t child_lidx = n[lidx].left_idx;
        uint32_t child_ridx = n[lidx].right_idx;
        printf("Min comparision ERROR; index: %u\n\
              \tparent: %u\n\
              \tflag:   %d\n\
              \tleft:   %u\n\
              \tmin:     %u, %lf\n\
              \tright:  %u\n\
              \tmin:     %u, %lf\n\
              \tnumPts: %u\n\
              \tmin:    %u\n\
              \tobject: %u\n",
              index,
              n[lidx].parent_idx, 
              flags[parent],
              child_lidx, 
              n[child_lidx].min_idx, 
              p[n[child_lidx].min_idx].z, 
              child_ridx, 
              n[child_ridx].min_idx, 
              p[n[child_ridx].min_idx].z, 
              n[lidx].numPts, 
              n[lidx].min_idx, 
              n[lidx].object_idx);
      }
#endif
      
      if(p[lmin].z - p[rmin].z < TOL)
        n[parent].min_idx = lmin;
      else
        n[parent].min_idx = rmin;

      // count points
      n[parent].numPts = n[lidx].numPts + n[ridx].numPts;
      // merge aabbs
      bb[parent] = merge(lbox, rbox);

      // look the next parent...
      parent = n[parent].parent_idx;
    }

    return;

  }

  return;
}

__global__
void create_aabb_size2( node_t* n, aabb_t* bb, int* flags, const point_t* p, 
  const size_t num_internal_nodes, const size_t num_nodes)
{
  uint32_t index = num_internal_nodes + threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i < num_nodes; i += stride)
  {
    uint32_t parent = n[i].parent_idx;

    while(parent != 0xFFFFFFFF && atomicCAS(flags + parent, 0, 1)) 
    {
      // here, the flag has already been 1. it means that this
      // thread is the 2nd thread. merge AABB of both childlen.

      const uint32_t lidx = n[parent].left_idx;
      const uint32_t ridx = n[parent].right_idx;
      const aabb_t lbox = bb[lidx];
      const aabb_t rbox = bb[ridx];
      const uint32_t lmin = n[lidx].min_idx;
      const uint32_t rmin = n[ridx].min_idx;
      // compare mins
      if(lmin == 0xFFFFFFFFu){
        uint32_t child_lidx = n[lidx].left_idx;
        uint32_t child_ridx = n[lidx].right_idx;
        printf("Min comparision ERROR\n\
              \tparent: %u\n\
              \tleft:   %u\n\
              \tmin:     %u, %lf\n\
              \tright:  %u\n\
              \tmin:     %u, %lf\n\
              \tnumPts: %u\n\
              \tmin:    %u\n\
              \tobject: %u\n", 
              n[lidx].parent_idx, 
              child_lidx, 
              n[child_lidx].min_idx, 
              p[n[child_lidx].min_idx].z, 
              child_ridx, 
              n[child_ridx].min_idx, 
              p[n[child_ridx].min_idx].z, 
              n[lidx].numPts, 
              n[lidx].min_idx, 
              n[lidx].object_idx);
      }
      
      if(p[lmin].z - p[rmin].z < TOL)
        n[parent].min_idx = lmin;
      else
        n[parent].min_idx = rmin;

      // count points
      n[parent].numPts = n[lidx].numPts + n[ridx].numPts;
      // merge aabbs
      bb[parent] = merge(lbox, rbox);

      // look the next parent...
      parent = n[parent].parent_idx;
    }

  }

  return;
}

// void create_aabb_cpu( const node_t* n, aabb_t* bb, const uint32_t current)
// {
//   if(n[current].object_idx != 0xFFFFFFFF)
//   {
//     return;
//   }
//   else
//   {
//     const uint32_t lidx = n[current].left_idx;
//     const uint32_t ridx = n[current].right_idx;
//     tbb::parallel_invoke([&]() { create_aabb_cpu( n, bb, lidx);}, [&]() { create_aabb_cpu( n, bb, ridx);} );
//     // create_aabb_cpu( n, bb, lidx);
//     // create_aabb_cpu( n, bb, ridx);
//     const aabb_t lbox = bb[lidx];
//     const aabb_t rbox = bb[ridx];
//     bb[current] = merge(lbox, rbox);
//   }
//   return;
// }


__global__
void init_struct(node_t* n, aabb_t* bb, const size_t NN, const node_t default_node, const aabb_t default_aabb)
{
  uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<NN; i += stride)
  {
    n[i] = default_node;
    bb[i] = default_aabb;
  }
  return;
}

// template<typename T, typename L, typename def>
// __global__
// void init_struct(T* pointer, const L length, const def default_item)
// {
//   uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
//   uint32_t stride = blockDim.x * gridDim.x;

//   for(uint32_t i = index; i<length; i += stride)
//   {
//     pointer[i] = default_item;
//   }
//   return;
// }



int main(int argc, char* argv[])
{
  std::string inputTXT = (argc > 1)? argv[1] : "data/INAER_2011_Alcoy.xyz";
  std::size_t N;
  aabb_t BBox;

  // node_t default_node;
  // default_node.parent_idx = 0xFFFFFFFF;
  // default_node.left_idx   = 0xFFFFFFFF;
  // default_node.right_idx  = 0xFFFFFFFF;
  // default_node.numPts     = 0xFFFFFFFF;
  // default_node.min_idx    = 0xFFFFFFFF;
  // default_node.object_idx = 0xFFFFFFFF;

  aabb_t default_aabb;
  default_aabb.upper.x = -inf; default_aabb.lower.x = inf;
  default_aabb.upper.y = -inf; default_aabb.lower.y = inf;
  // default_aabb.upper.z = -inf; default_aabb.lower.z = inf;

  // const std::size_t N = (argc > 1)? atoi(argv[1]) : 10;

  std::mt19937 mt(123456789);
  // std::uniform_real_distribution<real_t> uni(0.0, 1.0);
  std::uniform_real_distribution<real_t> uni(0.0, 10.0);

  if( inputTXT.find("INAER_2011_Alcoy.xyz") != std::string::npos ){
    N = 2772832;
    BBox.lower.x   = 715244.96;
    BBox.lower.y   = 4286623.63;
    // BBox.lower.z   = 836.424;
    BBox.upper.x   = 716057.75;
    BBox.upper.y   = 4287447.70;
    // BBox.upper.z   = 976.790;

  } else if( inputTXT.find("INAER_2011_Alcoy_Core.xyz") != std::string::npos ){ // Alcoy
    N = 20380212;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    // BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    // BBox.upper.z = 991.516;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX2.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*2;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    // BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    // BBox.upper.z = 991.516;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX4.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*4;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    // BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    // BBox.upper.z = 991.516;

  } else if( inputTXT.find("INAER_2011_Alcoy_CoreX6.xyz") != std::string::npos ){ // Alcoy
    N = 20380212*6;
    BBox.lower.x = 714947.98;
    BBox.lower.y = 4286501.93;
    // BBox.lower.z = 830.381;
    BBox.upper.x = 716361.06;
    BBox.upper.y = 4288406.23;
    // BBox.upper.z = 991.516;

  } else if( inputTXT.find("BABCOCK_2017_Arzua_3B.xyz") != std::string::npos ){ //Arzua
    N = 40706503;
    BBox.lower.x = 568000.00;
    BBox.lower.y = 4752320.00;
    // BBox.lower.z = 331.620;
    BBox.upper.x = 568999.99;
    BBox.upper.y = 4753319.99;
    // BBox.upper.z = 495.630;

  } else if( inputTXT.find("V21_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion forestal
    N = 42384876;
    BBox.lower.x = 526964.093;
    BBox.lower.y = 4742610.292;
    // BBox.lower.z = 38.656;
    BBox.upper.x = 527664.647;
    BBox.upper.y = 4743115.738;
    // BBox.upper.z = 112.269;

  } else if( inputTXT.find("V19_group1_densified_point_cloud.xyz") != std::string::npos ){ //Brion urban
    N = 48024480;
    BBox.lower.x = 526955.908;
    BBox.lower.y = 4742586.025;
    // BBox.lower.z = 38.150;
    BBox.upper.x = 527686.445;
    BBox.upper.y = 4743124.373;
    // BBox.upper.z = 119.833;

  } else {
    printf("No header data!\n");
    exit(-1);
  }

  N = (argc > 5)? static_cast<size_t>(atoi(argv[5])) : N; 

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, 8);

  real_v diffBox;
  diffBox.x = BBox.upper.x - BBox.lower.x;
  diffBox.y = BBox.upper.y - BBox.lower.y;
  // diffBox.z = BBox.upper.z - BBox.lower.z;


  std::vector<point_t> ps(N);

  if(read_pointsC(inputTXT, ps) < 0){
      printf("Unable to read file!\n");
      exit(-1);
  }

  point_t*   point_cloud = NULL;
  point_t*   point_cloud_out = NULL;
  
  uint32_t*  morton = NULL;
  // uint64_t*  morton64 = NULL;
  uint32_t*  morton_out = NULL;
  uint32_t*  indices = NULL;
  uint32_t*  indices_out = NULL;
  void*      d_temp_storage = NULL;

  // node_t* node_list = NULL;
  // node_t* aux_node_list = NULL;
  // aabb_t* aabb_list = NULL;
  // aabb_t* aux_aabb_list = NULL;

  // int* flags;
  int* flagsS;

  uint64_t* morton64S = NULL;
  node_t* node_listS = NULL;
  aabb_t* aabb_listS = NULL;

  num_objects = N;
  // const size_t num_internal_nodes = num_objects - 1;
  // const size_t num_nodes = 2*num_objects - 1; /*Numero de Nodos*/ 

  int deviceId;
  int numberOfSMs;

  const size_t size_points = num_objects * sizeof(point_t);
  const size_t size_morton = num_objects * sizeof(uint32_t);
  // const size_t size_morton64 = num_objects * sizeof(uint64_t);
  size_t temp_storage_bytes = 0;
  
  // const size_t size_nodes = num_nodes * sizeof(node_t);
  // const size_t size_aabbs = num_nodes * sizeof(aabb_t);

  cudaError_t mortonError;
  cudaError_t asyncErr;

  cudaGetDevice(&deviceId);
  cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);
  printf("Device ID: %d\tNumber of SMs: %d\n", deviceId, numberOfSMs);

  const int chunkDim = (argc > 2)? atoi(argv[2]) : 16;
  r_num_objects = (uint32_t)((num_objects-1)/chunkDim) + 1;
  r_num_internal_nodes = r_num_objects - 1;
  r_num_nodes = 2*r_num_objects -1;

  const uint32_t threadsPerBlock = (argc > 3)? static_cast<uint32_t>(atoi(argv[3])) : 256;
  const uint32_t numberOfBlocks = (argc > 4)? static_cast<uint32_t>(atoi(argv[4]))*numberOfSMs : 32*numberOfSMs;
  // size_t threadsPerBlock = 1;
  // size_t numberOfBlocks = 1 * numberOfSMs;
  const uint32_t nBlocks_node_func = (uint32_t)((r_num_internal_nodes-1)/threadsPerBlock) + 1;
  const uint32_t nBlocks_aabb_func = (uint32_t)((r_num_objects-1)/threadsPerBlock) + 1;

  std::cout << inputTXT << "; " << N << " points; chunkDim: " << chunkDim;
  std::cout << "; threadsPerBlock: " << threadsPerBlock << "; numBlocks: " << numberOfBlocks << std::endl;
  std::cout << r_num_objects << " leaf nodes\n";
  std::cout << r_num_internal_nodes << " internal nodes\n";

  cudaMallocManaged(&point_cloud, size_points);

  cudaMallocManaged(&morton, size_morton);
  // cudaMallocManaged(&morton64, size_morton64);
  cudaMallocManaged(&morton_out, size_morton);
  cudaMallocManaged(&indices, size_morton);
  cudaMallocManaged(&indices_out, size_morton);

  cudaMallocManaged(&morton64S, r_num_objects*sizeof(uint64_t));
  cudaMallocManaged(&node_listS, r_num_nodes*sizeof(node_t));
  cudaMallocManaged(&aabb_listS, r_num_nodes*sizeof(aabb_t));

  // cudaMallocManaged(&flags, N*sizeof(int));
  cudaMallocManaged(&flagsS, r_num_internal_nodes*sizeof(int));

  // cudaMemPrefetchAsync(flagsS, r_num_internal_nodes*sizeof(int), cudaCpuDeviceId);

  // for(int i=0; i<r_num_internal_nodes; i++){
  //   if(flagsS[i]!=0)
  //     printf("Error flags\n");
  //   flagsS[i] = 0;
  // }

  cudaMemPrefetchAsync(node_listS, r_num_nodes*sizeof(node_t), deviceId);
  cudaMemPrefetchAsync(aabb_listS, r_num_nodes*sizeof(aabb_t), deviceId);

  // cudaMemPrefetchAsync(node_list, size_nodes, deviceId);
  // cudaMemPrefetchAsync(aabb_list, size_aabbs, deviceId);
  cudaMemPrefetchAsync(morton, size_morton, deviceId);
  cudaMemPrefetchAsync(indices, size_morton, deviceId);
  cudaMemPrefetchAsync(morton_out, size_morton, deviceId);
  cudaMemPrefetchAsync(indices_out, size_morton, deviceId);
  // cudaMemPrefetchAsync(morton64, size_morton64, deviceId);
  cudaMemPrefetchAsync(morton64S, r_num_objects*sizeof(uint64_t), deviceId);
  // cudaMemPrefetchAsync(flags, num_objects*sizeof(int), deviceId);
  cudaMemPrefetchAsync(flagsS, r_num_internal_nodes*sizeof(int), deviceId);

  cudaMemPrefetchAsync(point_cloud, size_points, cudaCpuDeviceId);

#ifdef DEBUG
  int i = 0;
  for(auto& p : ps){
    point_cloud[i].x = (int)uni(mt);
    point_cloud[i].y = (int)uni(mt);
    point_cloud[i].z = (int)uni(mt);
    i++;
  }
  point_cloud[9].z = -1.0;
#else
  int i = 0;
  for(auto& p : ps){
    point_cloud[i].x = p.x;
    point_cloud[i].y = p.y;
    point_cloud[i].z = p.z;
    i++;
  }
#endif
  
  cudaMemPrefetchAsync(point_cloud, size_points, deviceId);
  // cudaMemPrefetchAsync(point_cloud_out, size_points, deviceId);

  // cudaMemPrefetchAsync(node_list, size_nodes, deviceId);
  // cudaMemPrefetchAsync(aabb_list, size_aabbs, deviceId);

  // init_struct<<<numberOfBlocks, threadsPerBlock>>>( node_list, aabb_list, num_nodes, default_node, default_aabb);
  // init_struct<<<numberOfBlocks, threadsPerBlock>>>( aabb_list, num_internal_nodes, default_aabb);
  // init_struct<<<numberOfBlocks, threadsPerBlock>>>( node_list, num_internal_nodes, default_node);

  init_struct<<<numberOfBlocks, threadsPerBlock>>>( node_listS, aabb_listS, r_num_nodes, default_node, default_aabb);

  mortonError = cudaGetLastError();
  if(mortonError != cudaSuccess) printf("Error INIT NODES y ABBBs: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error INIT NODES y ABBBs sync: %s\n", cudaGetErrorString(asyncErr));


  // cudaMemPrefetchAsync(point_cloud, size_points, deviceId);
  // cudaMemPrefetchAsync(morton, size_morton, deviceId);
  // cudaMemPrefetchAsync(indices, size_morton, deviceId);

  // double totaltime = 0;
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();

  get_morton_code<<<numberOfBlocks, threadsPerBlock>>>( point_cloud, BBox, diffBox, morton, indices, num_objects);

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error MORTON: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error MORTON sync: %s\n", cudaGetErrorString(asyncErr));

  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " MORTON time elapsed: " << mytime << " ms\n";
  double totaltime = mytime;

#ifdef DEBUG
  for(int i = 0; i<N; i++){
    std::cout << point_cloud[i].x << ",";
    std::cout << point_cloud[i].y << ",";
    std::cout << point_cloud[i].z << "  ";
  }
  std::cout << "\n";
  for(int i = 0; i<N; i++){
    std::cout << std::bitset<32>(morton[i]) << " ";
  }
  std::cout << "\n";
  for(int i = 0; i<N; i++){
    std::cout << indices[i] << " ";
  }
  std::cout << "\n";
#endif

  // // cudaMemPrefetchAsync(morton, size_morton, deviceId);
  // cudaMemPrefetchAsync(morton_out, size_morton, deviceId);
  // // cudaMemPrefetchAsync(indices, size_morton, deviceId);
  // cudaMemPrefetchAsync(indices_out, size_morton, deviceId);


  i_start = tempo_t::now();

  /* Determine temporary device storage requirements; segun las especificaciones
  si el puntero temporal apunta a NULL, se modifica "temp_storage_bytes" con el
  tamaño de memoria temporal requerida */

  cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
      morton, morton_out, indices, indices_out, N);
  /* Allocate temporary storage */
  cudaMallocManaged(&d_temp_storage, temp_storage_bytes);
  cudaMemPrefetchAsync(d_temp_storage, temp_storage_bytes, deviceId);
  
  /* solo ordeno los índices */
  cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
      morton, morton_out, indices, indices_out, N);

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error SORT: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error SORT sync: %s\n", cudaGetErrorString(asyncErr));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " SORT time elapsed: " << mytime << " ms\n";
  totaltime += mytime;

#ifdef DEBUG
  for(int i = 0; i<N; i++){
    std::cout << std::bitset<32>(morton_out[i]) << " ";
  }
  std::cout << "\n";
  for(int i = 0; i<N; i++){
    std::cout << indices_out[i] << " ";
  }
  std::cout << "\n";
  for(int i = 0; i<N; i++){
    std::cout << point_cloud[indices_out[i]].x << ",";
    std::cout << point_cloud[indices_out[i]].y << ",";
    std::cout << point_cloud[indices_out[i]].z << "  ";
  }
  std::cout << "\n";
#endif



//   // cudaMemPrefetchAsync(morton64, size_morton64, deviceId);

//   i_start = tempo_t::now();

//   // check_morton<<<numberOfBlocks, threadsPerBlock>>>( morton_out, indices_out, morton64, num_objects );

//   check_morton_size<<<numberOfBlocks, threadsPerBlock>>>( morton_out, indices_out, morton64S, r_num_objects, chunkDim);

//   // mortonError = cudaGetLastError();
//   // if(mortonError != cudaSuccess) printf("Error CHECK: %s\n", cudaGetErrorString(mortonError));

//   // asyncErr = cudaDeviceSynchronize();
//   // if(asyncErr != cudaSuccess) printf("Error CHECK sync: %s\n", cudaGetErrorString(asyncErr));

//   mytime = cast_t(tempo_t::now() - i_start).count();
//   std::cout << " CHECK time elapsed: " << mytime << " ms\n";
//   totaltime += mytime;

// #ifdef DEBUG
//   // for(int i = 0; i<N; i++){
//   //   std::cout << std::bitset<64>(morton64[i]) << " ";
//   // }
//   std::cout << "\n";
//   for(int i = 0; i<N; i++){
//     std::cout << indices_out[i] << " ";
//   }
//   std::cout << "\n";
//   for(int i = 0; i<r_num_objects; i++){
//     std::cout << std::bitset<64>(morton64S[i]) << " ";
//   }
//   std::cout << "\n";
// #endif


  // cudaMemPrefetchAsync(aux_node_list, num_objects*sizeof(node_t), deviceId);
  // cudaMemPrefetchAsync(aux_aabb_list, num_objects*sizeof(aabb_t), deviceId);

  i_start = tempo_t::now();

  // cudaMemsetAsync( node_list, 0xFFFFFFFF, size_nodes*0.5 , (cudaStream_t)1);

  // init_leafs<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, aux_aabb_list, point_cloud, indices_out, num_objects );

  // init_leafs_size<<<numberOfBlocks, threadsPerBlock>>>( &node_listS[r_num_internal_nodes], 
  //   &aabb_listS[r_num_internal_nodes], 
  //   point_cloud, 
  //   indices_out, 
  //   num_objects,
  //   r_num_objects,
  //   chunkDim,
  //   default_aabb );

  init_leafs_size2<<<numberOfBlocks, threadsPerBlock>>>( &node_listS[r_num_internal_nodes], 
    &aabb_listS[r_num_internal_nodes], 
    point_cloud, 
    indices_out,
    morton64S, 
    num_objects,
    r_num_objects,
    chunkDim,
    default_aabb,
    BBox,
    diffBox);
  
  // init_leafs_size<<<numberOfBlocks, threadsPerBlock>>>( &node_listS[r_num_internal_nodes], 
  //                                                       &aabb_listS[r_num_internal_nodes], 
  //                                                       point_cloud,
  //                                                       point_cloud_out, 
  //                                                       indices_out, 
  //                                                       num_objects,
  //                                                       r_num_objects,
  //                                                       chunkDim,
  //                                                       default_aabb );

// init_leafs<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, indices, N );

  // init_leafs2<<<numberOfBlocks, threadsPerBlock>>>( node_list, aabb_list, point_cloud, indices, num_internal_nodes, num_nodes );

  // init_indices<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, indices, N );

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error LEAFS: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error LEAFS sync: %s\n", cudaGetErrorString(asyncErr));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " LEAFS and AABBs time elapsed: " << mytime << " ms\n";
  totaltime += mytime;

#ifdef DEBUG
  for(int i = 0; i<r_num_nodes; i++) {
    std::cout << node_listS[i].parent_idx << ", " << node_listS[i].left_idx << ", ";
    std::cout << node_listS[i].right_idx << ", " << node_listS[i].object_idx << ", ";
    std::cout << node_listS[i].numPts << ", " << node_listS[i].min_idx << "\n";
  }
  std::cout << std::endl;

  for(int i = 0; i<r_num_nodes; i++){
    std::cout << aabb_listS[i].upper.x << "," << aabb_listS[i].upper.y << " ";
    std::cout << aabb_listS[i].lower.x << "," << aabb_listS[i].lower.y << "\n";
  }
  std::cout << std::endl;
#endif

  // cudaMemPrefetchAsync(node_list, size_nodes, deviceId);
  // cudaMemPrefetchAsync(morton64, size_morton64, deviceId);

  i_start = tempo_t::now();

  // init_nodes<<<numberOfBlocks, threadsPerBlock>>>( node_listS, morton64S, r_num_objects );
  init_nodes2<<<nBlocks_node_func, threadsPerBlock>>>( node_listS, morton64S, r_num_objects );

//   thrust::for_each(thrust::device,
//     thrust::make_counting_iterator<uint32_t>(0),
//     thrust::make_counting_iterator<uint32_t>(num_internal_nodes),
//     [node_list, morton64, num_objects] __device__ (const uint32_t idx)
//     {
//       // node_t aux;

//       const uint2 ij  = lbvh::detail::determine_range(morton64, num_objects, idx);
//       // const uint32_t gamma = lbvh::detail::find_split(morton64, num_objects, ij.x, ij.y);
//       uint32_t left_idx = lbvh::detail::find_split(morton64, num_objects, ij.x, ij.y);

// // #ifdef DEBUG
// //       printf("%u, %u left_idx: %d\n",ij.x, ij.y, left_idx);
// // #endif

//       // aux.object_idx = 0xFFFFFFFF;

//       // node_list[idx].left_idx  = gamma;
//       // node_list[idx].right_idx = gamma + 1;
//       uint32_t right_idx = left_idx + 1;

//       if(thrust::min(ij.x, ij.y) == left_idx)
//       {
//           left_idx += num_objects - 1;
//       }
//       if(thrust::max(ij.x, ij.y) == right_idx)
//       {
//           right_idx += num_objects - 1;
//       }
//       // node_list[idx] = aux;

//       node_list[idx].left_idx = left_idx;
//       node_list[idx].right_idx = right_idx;
//       node_list[left_idx].parent_idx  = idx;
//       node_list[right_idx].parent_idx = idx;
//       return;
//     });

  // thrust::for_each(thrust::device,
  //     thrust::make_counting_iterator<uint32_t>(0),
  //     thrust::make_counting_iterator<uint32_t>(r_num_internal_nodes),
  //     [node_listS, morton64S, r_num_objects] __device__ (const uint32_t idx)
  //     {
  //       // node_t aux;

  //       const uint2 ij  = lbvh::detail::determine_range(morton64S, r_num_objects, idx);
  //       // const uint32_t gamma = lbvh::detail::find_split(morton64S, r_num_objects, ij.x, ij.y);
  //       uint32_t left_idx = lbvh::detail::find_split(morton64S, r_num_objects, ij.x, ij.y);

  //       uint32_t right_idx = left_idx + 1;

  //       if(thrust::min(ij.x, ij.y) == left_idx)
  //       {
  //           left_idx += r_num_objects - 1;
  //       }
  //       if(thrust::max(ij.x, ij.y) == right_idx)
  //       {
  //           right_idx += r_num_objects - 1;
  //       }
  //       // node_listS[idx] = aux;

  //       node_listS[idx].left_idx = left_idx;
  //       node_listS[idx].right_idx = right_idx;
  //       node_listS[left_idx].parent_idx  = idx;
  //       node_listS[right_idx].parent_idx = idx;
  //       return;
  //     });

// mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error NODES: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error NODES sync: %s\n", cudaGetErrorString(asyncErr));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " INTERNAL NODES time elapsed: " << mytime << " ms\n";
  totaltime += mytime;

// #ifdef DEBUG
//   for(int i = 0; i<r_num_nodes; i++){
//     std::cout << node_listS[i].parent_idx << ", " << node_listS[i].left_idx << ", ";
//     std::cout << node_listS[i].right_idx << ", " << node_listS[i].object_idx << ", ";
//     std::cout << node_listS[i].numPts << ", " << node_listS[i].min_idx << "\n";
//   }
//   std::cout << std::endl;
//   for(int i = 0; i<r_num_nodes; i++){
//       std::cout << aabb_listS[i].upper.x << "," << aabb_listS[i].upper.y << " ";
//       std::cout << aabb_listS[i].lower.x << "," << aabb_listS[i].lower.y << "\n";
//     }
//     std::cout << std::endl;
// #endif
  
  // cudaMemPrefetchAsync(flags, num_objects*sizeof(int), deviceId);
  // cudaMemPrefetchAsync(aabb_list, size_aabbs, deviceId);

  i_start = tempo_t::now();

  // uint32_t nBlocks_aabb_func = (uint32_t)(N/threadsPerBlock) + 1;

  // create_aabb<<<nBlocks_aabb_func, threadsPerBlock>>>( node_list, aabb_list, flags, num_internal_nodes, num_nodes );
  // create_aabb2<<<numberOfBlocks, threadsPerBlock>>>( node_list, aabb_list, flags, num_internal_nodes, num_nodes );
  // create_aabb_cpu(node_list, aabb_list, 0);

  create_aabb_size<<<nBlocks_aabb_func, threadsPerBlock>>>( node_listS, 
    aabb_listS, 
    flagsS,
    point_cloud, 
    r_num_internal_nodes, 
    r_num_nodes );
  // create_aabb_size2<<<numberOfBlocks, threadsPerBlock>>>( node_listS, 
  //   aabb_listS, 
  //   flagsS,
  //   point_cloud, 
  //   r_num_internal_nodes, 
  //   r_num_nodes );

  // thrust::for_each(thrust::device,
  //     thrust::make_counting_iterator<uint32_t>(r_num_internal_nodes),
  //     thrust::make_counting_iterator<uint32_t>(r_num_nodes),
  //     [node_listS, aabb_listS, point_cloud, flagsS] __device__ (const uint32_t idx)
  //     {
  //         uint32_t parent = node_listS[idx].parent_idx;
  //         while(parent != 0xFFFFFFFF) // means idx == 0
  //         {
  //           const int old = atomicCAS(flagsS + parent, 0, 1);
  //           if(old == 0)
  //           {
  //               // this is the first thread entered here.
  //               // wait the other thread from the other child node.
  //               return;
  //           }
  //           assert(old == 1);
  //           // here, the flag has already been 1. it means that this
  //           // thread is the 2nd thread. merge AABB of both childlen.
      
  //           const uint32_t lidx = node_listS[parent].left_idx;
  //           const uint32_t ridx = node_listS[parent].right_idx;
  //           const aabb_t lbox = aabb_listS[lidx];
  //           const aabb_t rbox = aabb_listS[ridx];
  //           const uint32_t lmin = node_listS[lidx].min_idx;
  //           const uint32_t rmin = node_listS[ridx].min_idx;
  //           // compare mins
            
  //           if(point_cloud[lmin].z - point_cloud[rmin].z < TOL)
  //             node_listS[parent].min_idx = lmin;
  //           else
  //             node_listS[parent].min_idx = rmin;
      
  //           // count points
  //           node_listS[parent].numPts = node_listS[lidx].numPts + node_listS[ridx].numPts;
  //           // merge aabbs
  //           aabb_listS[parent] = merge(lbox, rbox);
      
  //           // look the next parent...
  //           parent = node_listS[parent].parent_idx;
  //         }
  //         return;
  //     });

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error AABBs: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error AABBs sync: %s\n", cudaGetErrorString(asyncErr));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " Create AABB time elapsed: " << mytime << " ms\n";
  totaltime += mytime;
  std::cout << "  CREATION takes: " << totaltime << " ms\n";

#ifdef DEBUG
  for(int i = 0; i<r_num_nodes; i++){
    std::cout << node_listS[i].parent_idx << ", " << node_listS[i].left_idx << ", ";
    std::cout << node_listS[i].right_idx << ", " << node_listS[i].object_idx << ", ";
    std::cout << node_listS[i].numPts << ", " << node_listS[i].min_idx << "(";
    if(node_listS[i].min_idx != 0xFFFFFFFF)
      std::cout << point_cloud[node_listS[i].min_idx].z;
    std::cout << ")\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<r_num_nodes; i++){
    std::cout << aabb_listS[i].upper.x << "," << aabb_listS[i].upper.y << " ";
    std::cout << aabb_listS[i].lower.x << "," << aabb_listS[i].lower.y << "\n";
  }
  std::cout << std::endl;
  for(int i=0; i<r_num_internal_nodes; i++){
    std::cout << flagsS[i] << " ";
  }
  std::cout << std::endl;
#endif


#ifdef CHECK
  std::cout << "\n\n";
  lbvh::bvh<real_t, real_v, aabb_getter> bvh(ps.begin(), ps.end(), true);
  std::cout << "\n\n";


  auto bvh_nodes = bvh.nodes_host();
  auto bvh_aabbs = bvh.aabbs_host();
  auto bvh_morton = bvh.morton_host();
  auto bvh_morton64 = bvh.morton64_host();
  auto bvh_indices = bvh.indices_host();

  int error=0;
  for(int i=0; i<num_nodes; i++){
      // std::cout << bvh_nodes[i].parent_idx << ", " << node_list[i].parent_idx << "\n";
      if(bvh_nodes[i].parent_idx != node_list[i].parent_idx ||
          bvh_nodes[i].left_idx != node_list[i].left_idx ||
          bvh_nodes[i].right_idx != node_list[i].right_idx ||
          bvh_nodes[i].object_idx != node_list[i].object_idx )
        error++;
  }

  if(error)
    std::cout << "There are DIFFERENCES in NODES: " << error << std::endl;
  else
    std::cout << "NODES are OK!\n";
  
  error = 0;
  for(int i=0; i<num_nodes; i++){
      if( fabs(bvh_aabbs[i].upper.x - aabb_list[i].upper.x) > TOL ||
          fabs(bvh_aabbs[i].upper.y - aabb_list[i].upper.y) > TOL ||
          fabs(bvh_aabbs[i].upper.z - aabb_list[i].upper.z) > TOL ||
          fabs(bvh_aabbs[i].lower.x - aabb_list[i].lower.x) > TOL ||
          fabs(bvh_aabbs[i].lower.y - aabb_list[i].lower.y) > TOL ||
          fabs(bvh_aabbs[i].lower.z - aabb_list[i].lower.z) > TOL )
        error++;
  }

  if(error)
    std::cout << "There are DIFFERENCES in AABBS: " << error << std::endl;
  else
    std::cout << "AABBs are OK!\n";


  error=0;
  for(int i=0; i<num_objects; i++){
      if(bvh_morton[i] != morton_out[i])
        error++;
  }

  if(error)
    std::cout << "There are DIFFERENCES in MORTON: " << error << std::endl;
  else
    std::cout << "MORTONs are OK!\n";

  error=0;
  for(int i=0; i<num_objects; i++){
      if(bvh_indices[i] != indices_out[i])
        error++;
  }

  if(error)
    std::cout << "There are DIFFERENCES in INDICES: " << error << std::endl;
  else
    std::cout << "INDICES are OK!\n";


  error=0;
  for(int i=0; i<num_objects; i++){
      if(bvh_morton64[i] != morton64[i])
        error++;
  }

  if(error)
    std::cout << "There are DIFFERENCES in MORTON64: " << error << std::endl;
  else
    std::cout << "MORTON64s are OK!\n";


  bvh.clear();

#endif

  cudaFree(morton);
  cudaFree(morton_out);
  cudaFree(indices);
  cudaFree(d_temp_storage);
  cudaFree(morton64S);
  cudaFree(flagsS);


  // aabb_t mybb;
  // // mybb.upper.x = BBox.upper.x-diffBox.x*0.5;
  // mybb.upper.x = BBox.upper.x;
  // mybb.upper.y = BBox.upper.y;
  // mybb.lower.x = BBox.lower.x;
  // mybb.lower.y = BBox.lower.y;
  // mybb.upper.x = 5;
  // mybb.upper.y = 10;
  // mybb.lower.x = 0;
  // mybb.lower.y = 0;

  /*ahora compruebo si está bien construido el árbol*/
  // uint32_t numPts = 0;
  // point_t min = findValidMin(node_listS, aabb_listS, point_cloud, indices_out, node_listS[0], mybb, numPts);
  // if(node_listS[0].numPts != numPts){
  //   std::cout << "Count: " << node_listS[0].numPts << " != " << numPts << "; ";
  //   std::cout << node_listS[0].numPts - numPts << " less points\n";
  // }
  // else
  //   std::cout << "Count: " << node_listS[0].numPts << " == " << numPts << std::endl;

  // std::cout << "The minimum is: " << min.z << std::endl;


  // cudaMemPrefetchAsync(node_listS, r_num_nodes*sizeof(node_t), cudaCpuDeviceId);
  // cudaMemPrefetchAsync(aabb_listS, r_num_nodes*sizeof(aabb_t), cudaCpuDeviceId);
  // cudaMemPrefetchAsync(indices_out, size_morton, cudaCpuDeviceId);
  // cudaMemPrefetchAsync(point_cloud, size_points, cudaCpuDeviceId);

  uint32_t Wsize = 10;
  uint32_t Bsize = 20;
  double Overlap = 0.8;

  uint32_t Ncells, Ngrid;
  uint32_t nRows, nCols, nColsg, nRowsg;
  double Width, High, Density;

  // double best_time = 111111.0;
  // double best_rate;

  uint32_t countMin = 0;

  // uint32_t searcher = 0;

  // uint32_t addMin = 0;

  Width = diffBox.x;
  High = diffBox.y;
  // Densidad en puntos/m^2
  Density = N/(Width*High);

  printf("CLOUD PARAMETERS:\n");
  printf("Número de puntos      %lu\n",N);
  printf("Ancho:  %.2lf\n",Width);
  printf("Alto:  %.2lf\n",High);
  printf("Densidad:  %.3lf\n",Density);

  printf("\nTamaño de ventana     %u\n", Wsize);
  printf("Tamaño de rejilla     %u\nSolape                %.2f\n", Bsize,Overlap);

  // El numero minimo sera la mitad del numero de puntos medio por celda
  uint32_t minNumPoints = (uint32_t)(0.5*Density*Wsize*Wsize);
  printf("Minimo numero de puntos por celda:   %u\n", minNumPoints);

  double Displace = round2d(Wsize*(1-Overlap));
  printf("Displacement   %.2f\n", Displace);

  // Stage 1 parameters
  printf("\nVENTANA:\n");
  if(Overlap > 0.0) {
    nCols=(int)(round((Width+2*Wsize*Overlap)/Displace))-1;
    nRows=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
  } else {
    nCols=(int)floor(Width/Wsize)+1;
    nRows=(int)floor(High/Wsize)+1;
  }
  printf("Celdas por columa: %d\n",nRows);
  printf("Celdas por fila:   %d\n",nCols);
  Ncells = nCols*nRows;

  // Stage 3 parameters
  printf("\nMALLA:\n");
  nColsg=(int)floor(Width/Bsize)+1;
  nRowsg=(int)floor(High/Bsize)+1;
  printf("Dimensiones de la malla %dx%d\n\n", nRowsg, nColsg);
  Ngrid = nColsg*nRowsg;

  // double initX = BBox.lower.x - Wsize/2 + Displace;
  // double initY = BBox.lower.y - Wsize/2 + Displace;

  /* Este es el bloque que donde creo un árbol con la insercción de puntos normal
  y después copi el arbol utilizando tree2_t */
  // Vector2D min = {BBox.lower.x,BBox.lower.y};
  // Vector2D max = {BBox.upper.x,BBox.upper.y};
  // float maxRadius = 0.0;
  // Vector2D radius = getRadius(min, max, &maxRadius);
  // Vector2D center = getCenter(min, radius);
  // Btree cpu_btree = new Btree_t( NULL, center, maxRadius);
  // uint maxSize = 16;
  // for(int i = 0; i < N; i++){
  //   // insertPoint(&point_cloud[i], cpu_qtree, minRadius);
  //   insertPoint2(&point_cloud[i], cpu_btree, maxSize);
  // }
  // array_all_points = (point_t*)(std::malloc(size_points));
  // array_nodes2 = (tree2)(std::malloc(cpu_tree_nodes*sizeof(tree2_t)));
  // tree2 cpu_tree2 = launchCopy4(cpu_btree);
  // printf("El arbol original tiene %u nodos y yo creo %u alineados\n", cpu_tree_nodes, copied_nodes);


  /*Aquí hago el recorrido sobre el árbol GPU utilizando la CPU con el mismo método
  que utilizaba en el código del TFM*/
  point_t* minVector = static_cast<point_t*>(malloc( Ncells*sizeof(point_t) ));
  memset(minVector, 0, 4*Ncells*sizeof(double));
  countMin = 0;

  for(int i=0; i<4; i++){
    i_start = tempo_t::now();

    stage1tbbRem(node_listS, aabb_listS, point_cloud, indices_out, 
        Wsize, Overlap, nCols, nRows, minNumPoints, minVector, node_listS[0], BBox);
  
    mytime = cast_t(tempo_t::now() - i_start).count();
    std::cout << " Stage1 BARE TREE time elapsed: " << mytime << " ms\n";
  
    countMin=0;
    for(int i=0; i<Ncells; i++){
      if(minVector[i].x != 0)
        countMin++;
    }
    printf("Numero de minimos: %u\n", countMin);
  
    memset(minVector, 0, 4*Ncells*sizeof(double));
  }    
  free(minVector);
  


  printf("\n");
  /*Este es donde ya lanzo un kernel CUDA de búsqueda*/
  // aabb_t* list = NULL;
  uint32_t* count = NULL;
  // cudaMallocManaged(&list, Ncells*sizeof(aabb_t));
  cudaMallocManaged(&count, Ncells*sizeof(uint32_t));
  for(int i=0; i<Ncells; i++){
    count[i] = 0u;
  }

  for(int i=0; i<4; i++){
    cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), deviceId);

    i_start = tempo_t::now();

    stage1tbb(node_listS, aabb_listS, point_cloud, indices_out, count,
        Wsize, Overlap, nCols, nRows, minNumPoints, BBox,
        threadsPerBlock, r_num_internal_nodes);
  
    mytime = cast_t(tempo_t::now() - i_start).count();
    std::cout << " Stage1 KERNEL CUDA time elapsed: " << mytime << " ms\n";
  
    countMin=0;
    for(int i=0; i<Ncells; i++){
      if(count[i] != 0)
        countMin++;
      count[i] = 0u;
    }
    printf("Numero de minimos: %u\n", countMin);
  
    // cudaMemset((void*)countMin, 0u, Ncells*sizeof(uint32_t));
  }
  // cudaFree(list);
  // cudaFree(count);




  /*En este lo que hago es alinear el vector de nodos, AABBs y puntos para GPU*/
  cudaMallocManaged(&array_all_points, num_objects*sizeof(point_t));
  cudaMallocManaged(&array_nodes_gpu, r_num_nodes*sizeof(node_t));
  cudaMallocManaged(&array_aabbs_gpu, r_num_nodes*sizeof(aabb_t));

  cudaMemPrefetchAsync(array_all_points, num_objects*sizeof(point_t), cudaCpuDeviceId);
  cudaMemPrefetchAsync(array_nodes_gpu, r_num_nodes*sizeof(node_t), cudaCpuDeviceId);
  cudaMemPrefetchAsync(array_aabbs_gpu, r_num_nodes*sizeof(aabb_t), cudaCpuDeviceId);

  for(int i=0; i<r_num_nodes; i++){
    array_nodes_gpu[i] = default_node;
    array_aabbs_gpu[i] = default_aabb;
  }

  launchCopyGPU(node_listS, aabb_listS, point_cloud, indices_out);
  printf("\nEl arbol original tiene %u nodos y yo creo %u alineados\n", r_num_nodes, gpu_copied_nodes);
  printf("Root: %u puntos; Total copiado: %u puntos\n", array_nodes_gpu[0].numPts, num_pts_copiados);

  // uint32_t* count = NULL;
  // cudaMallocManaged(&count, Ncells*sizeof(uint32_t));
  // for(int i=0; i<Ncells; i++){
  //   count[i] = 0u;
  // }

  cudaMemPrefetchAsync(array_all_points, num_objects*sizeof(point_t), deviceId);
  cudaMemPrefetchAsync(array_nodes_gpu, r_num_nodes*sizeof(node_t), deviceId);
  cudaMemPrefetchAsync(array_aabbs_gpu, r_num_nodes*sizeof(aabb_t), deviceId);


  for(int i=0; i<4; i++){
    cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), deviceId);

    cudaError_t check = cudaGetLastError();
    if(check != cudaSuccess) printf("Error PREFETCH: %s\n", cudaGetErrorString(asyncErr));
    asyncErr = cudaDeviceSynchronize();
    if(asyncErr != cudaSuccess) printf("Error PREFETCH sync: %s\n", cudaGetErrorString(asyncErr));
  
    i_start = tempo_t::now();

    stage1tbb( array_nodes_gpu, array_aabbs_gpu, array_all_points, count,
      Wsize, Overlap, nCols, nRows, minNumPoints, BBox, threadsPerBlock );
    // stage1tbb(node_listS, aabb_listS, point_cloud, indices_out, count,
    //   Wsize, Overlap, nCols, nRows, minNumPoints, BBox, threadsPerBlock, r_num_internal_nodes);

    mytime = cast_t(tempo_t::now() - i_start).count();
    std::cout << " Stage1 KERNEL CUDA ALIGNED time elapsed: " << mytime << " ms\n";
  
    countMin=0;
    for(int i=0; i<Ncells; i++){
      if(count[i] != 0)
        countMin++;
      count[i] = 0u;
    }
    printf("Numero de minimos: %u\n", countMin);
  
    // cudaMemset((void*)countMin, 0u, Ncells*sizeof(uint32_t));
  }
  // cudaFree(list);
  cudaFree(count);




  // /*Aquí hago una copia alineada en CPU del árbol creado en GPU
  // para hacer la búsqueda SOLO CPU*/
  // array_all_points = (point_t*)(std::malloc(size_points));
  // array_nodes = (tree)(std::malloc(r_num_nodes*sizeof(tree_t)));
  // // array_nodes2 = (tree2)(std::malloc(r_num_nodes*sizeof(tree2_t)));
  // tree cpu_tree = launchCopy2(node_listS, aabb_listS, point_cloud, indices_out);
  // // tree2 cpu_tree2 = launchCopy3(node_listS, aabb_listS, point_cloud, indices_out);
  // printf("El arbol original tiene %u nodos y yo creo %u alineados\n", r_num_nodes, copied_nodes);
  // printf("Root: %u puntos; Total: %u puntos\n", cpu_tree->numPts, num_pts_copiados);
      
  // for(int i=0; i<4; i++){
  //   i_start = tempo_t::now();

  //   stage1tbbRem(Wsize, Overlap, nCols, nRows, minNumPoints, minVector, cpu_tree, BBox);

  //   mytime = cast_t(tempo_t::now() - i_start).count();
  //   std::cout << " Stage1 ALIGNED TREE time elapsed: " << mytime << " ms\n";
        
  //   countMin=0;
  //   for(int i=0; i<Ncells; i++){
  //     if(minVector[i].x != 0)
  //       countMin++;
  //   }
  //   printf("Numero de minimos: %u\n", countMin);

  //   memset(minVector, 0, 4*Ncells*sizeof(double));
  // }





  // for(int i=0; i<4; i++){
  //   i_start = tempo_t::now();

  //   stage1tbbRem2(Wsize, Overlap, nCols, nRows, minNumPoints, minVector, cpu_tree2, BBox.lower);

  //   mytime = cast_t(tempo_t::now() - i_start).count();
  //   std::cout << " Stage1 time elapsed: " << mytime << " ms\n";
        
  //   countMin=0;
  //   for(int i=0; i<Ncells; i++){
  //     if(minVector[i].x != 0)
  //       countMin++;
  //   }
  //   printf("Numero de minimos: %u\n", countMin);

  //   memset(minVector, 0, 4*Ncells*sizeof(double));
  // }


  // std::vector<int> list(24,0);
  // checkStruct(array_nodes_gpu, array_nodes_gpu[0], 0, list);

  // printf("Levels CUDA: %zu\n", list.size());
  // for (uint32_t level = 0; level < list.size(); ++level)
  //   printf("  level %u : %u nodes\n", level, list[level] );

  // memset(list.data(), 0, 24*sizeof(int));

  // checkStruct(node_listS, node_listS[0], 0, list);

  // printf("Levels CUDA: %zu\n", list.size());
  // for (uint32_t level = 0; level < list.size(); ++level)
  //   printf("  level %u : %u nodes\n", level, list[level] );

  // checkStruct(cpu_btree, 0, list);

  // printf("Levels CPU: %zu\n", list.size());
  // for (uint32_t level = 0; level < list.size(); ++level)
  //   printf("  level %u : %u nodes\n", level, list[level] );
    
  // memset(list.data(), 0, 24*sizeof(int));

  // checkStruct(cpu_tree, 0, list);

  // printf("Levels COPY: %zu\n", list.size());
  // for (uint32_t level = 0; level < list.size(); ++level)
  //   printf("  level %u : %u nodes\n", level, list[level] );


  // std::vector<std::pair<Vector2D,float>> areas;
  // checkArea(cpu_tree2, areas);
  // if(save_areas("areas.xyz", inputTXT, areas) < 0){
  //   printf("Unable to create results file!\n");
  // }
        
        
  /*ejecucion del código del TFM*/

  // Vector2D min = {BBox.lower.x,BBox.lower.y};
  // Vector2D max = {BBox.upper.x,BBox.upper.y};
  // float maxRadius = 0.0;

  // Vector2D radius = getRadius(min, max, &maxRadius);
  // Vector2D center = getCenter(min, radius);

  // Lpoint* point_cloud2 = static_cast<Lpoint*>(malloc(N * sizeof(Lpoint)));

  // if(read_pointsC(inputTXT, point_cloud2, N) < 0){
  //   printf("Unable to read file!\n");
  //   // exit(-1);
  // }

  // printf("INSERTING POINTS...\n");
  // i_start = tempo_t::now();

  // int max_level = 3;
  // int maxSize = 32;

  // Qtree cpu_qtree = parallel_qtree_pf2(max_level, center, maxRadius, maxSize, point_cloud2, N);

  // double cpu_tree_time = cast_t(tempo_t::now() - i_start).count()/1e3;
  // std::cout << "  INSERT CPU time elapsed: " << cpu_tree_time << " s\n";

  // array_pointers = static_cast<QtreeG4>(malloc( cpu_tree_nodes * sizeof(QtreeG4_t)));

  // array_all_points = static_cast<Lpoint*>(malloc( N * sizeof(Lpoint)));

  // i_start = tempo_t::now();

  // QtreeG4 aligned_qtree = launchCopy(cpu_qtree);

  // double new_tree_time = cast_t(tempo_t::now() - i_start).count()/1e3;
  // std::cout << "  COPY TREE time elapsed: " << new_tree_time << " s\n\n";

  // int* minIDs = static_cast<int*>(malloc( Ncells*sizeof(int) ));
  // memset(minIDs, -1, Ncells*sizeof(int));

  // i_start = tempo_t::now();

  // stage1tbbRem(Wsize, Overlap, nCols, nRows, minNumPoints, minIDs, aligned_qtree, min);

  // mytime = cast_t(tempo_t::now() - i_start).count();
  // std::cout << " Stage1 time elapsed: " << mytime << " ms\n";

  // countMin=0;
  // for(int i=0; i<Ncells; i++){
  //   if(minIDs[i] != -1)
  //     countMin++;
  // }
  // printf("Numero de minimos: %u\n", countMin);



  cudaFree(point_cloud);
  cudaFree(point_cloud_out);
  // cudaFree(morton);
  // cudaFree(morton64);
  // cudaFree(morton_out);
  // cudaFree(indices);
  cudaFree(indices_out);
  // cudaFree(d_temp_storage);

  // cudaFree(node_list);
  // aux_node_list = NULL;
  // cudaFree(aabb_list);
  // aux_aabb_list = NULL;

  // cudaFree(flags);
  // cudaFree(flagsS);

  // cudaFree(morton64S);
  cudaFree(node_listS);
  cudaFree(aabb_listS);

  // free(array_all_points);
  // free(array_nodes);
  // free(array_nodes2);

  // free(array_pointers);
  // free(minIDs);
  // deleteQtree(cpu_qtree);
  // delete(cpu_qtree);

  // free(point_cloud2);

  // deleteBtree(cpu_btree);
  // delete(cpu_btree);

  cudaFree(array_all_points);
  cudaFree(array_nodes_gpu);
  cudaFree(array_aabbs_gpu);


  return 0;
}
