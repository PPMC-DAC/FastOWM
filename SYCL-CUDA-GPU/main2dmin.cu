#include "lbvh.cuh"
#include <nppdefs.h>
#include <random>
#include <vector>
#include <chrono>
#include <fstream>
#include <thrust/random.h>
// #include <tbb/parallel_invoke.h>
#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include "tbb/global_control.h"
#include "tbb/parallel_for.h"
// #include "tbb/blocked_range2d.h"

// #include <cub/cub.cuh>
// #include <cub/device/device_radix_sort.cuh>
#include <bitset>
// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>


using real_t = double;
using real_v = double2;
using point_t = double4;
// using node_t = lbvh::detail::node;

const auto inf = std::numeric_limits<real_t>::infinity();

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

struct aabb_t
{
    real_v upper;
    real_v lower;
};

#define TOL 1e-10
#define FW_S32_MIN  (~0x7FFFFFFF)
#define NUM_PROCS 8

void handle_error(const char* msg) {
  perror(msg); 
  exit(255);
}

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

void map_file_atof_tbb_th(std::string fname, uint32_t num_objects, point_t* cloud)
{
    size_t length = 0u;
    int fd = open(fname.c_str(), O_RDONLY);
    if (fd < 0)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) < 0)
        handle_error("fstat");

    length = sb.st_size;

    // creates a new mapping in the virtual address space of the calling process
    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");

    // fd can be closed immediately without invalidating the mapping
    if(close(fd) < 0)
        handle_error("close");

    int stride = int(num_objects/NUM_PROCS) + 1;

    tbb::parallel_for( 0, NUM_PROCS, 1,
        [&](int id){

            int myline = stride * id;
            int fin = (stride*(id+1) <= num_objects)? stride*(id+1) : num_objects;
            auto mylist = addr;
            uint32_t m_numLines = 0;
            while(m_numLines < fin) {
                if(m_numLines == myline){
                    for(int i = myline; i<fin; i++){
                        cloud[i].x = atof(mylist);
                        cloud[i].y = atof(mylist+10);
                        cloud[i].z = atof(mylist+22);

                        mylist = static_cast<const char*>(memchr(mylist, '\n', 64)) + 1;
                    }
                    m_numLines += stride; //fin
                } 
                else {
                    mylist = static_cast<const char*>(memchr(mylist, '\n', 64)) + 1;
                    m_numLines++;
                }
                
            }
    });


    // unmap
    if(munmap((void*)addr, length) < 0)
        handle_error("munmap");

    return;
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
inline uint32_t morton2D(point_t xyz , const real_t resolution = 65536.0 ) noexcept
{
  xyz.x = fmin(xyz.x * resolution, resolution - 1.0);
  xyz.y = fmin(xyz.y * resolution, resolution - 1.0);

  const uint32_t xx = SeparateBy1(static_cast<uint32_t>(xyz.x));
  const uint32_t yy = SeparateBy1(static_cast<uint32_t>(xyz.y));
  return xx * 2 + yy;
}

__device__
inline uint32_t morton2D(const real_v& xy /*, const real_t resolution = 65536.0 */) noexcept
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

  /*el problema de hacer esto es que no estoy almacenando
  el aabb correcto que represental al conjunto de puntos,
  sino que estoy dando solo la posición del primerpunto
  de cada conjunto*/

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
  const int num_objects, const int r_objects, const int chunkDim, const aabb_t default_aabb,
  const aabb_t BBox, const real_v diffBox)
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  const uint32_t stride = blockDim.x * gridDim.x;
  point_t auxp;

  for(uint32_t i = index; i<r_objects; i += stride)
  {
    uint32_t start = i*chunkDim;
    uint32_t end = i*chunkDim+chunkDim;
    aabb_t auxbb = default_aabb;

    real_t min =  99999.0;
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

    real_v c = centroid(auxbb);

    c.x -= BBox.lower.x;
    c.y -= BBox.lower.y;
    c.x /= diffBox.x;
    c.y /= diffBox.y;

    uint64_t aux64 = morton2D(c);
    aux64 <<= 32;
    aux64 |= i;
    m64[i] = aux64;

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

__device__
inline void swap(uint32_t& left, uint32_t& right)
{
  uint32_t aux = left;
  left = right;
  right = aux;
}

__device__ __forceinline__
uint2 determine_range2(const uint64_t* node_code,
        const uint32_t num_leaves, uint32_t idx)
{
  if(idx==0) return make_uint2(0, num_leaves-1);

  // Choose direction.
  unsigned int code = node_code[idx];
  int prefix_prev = ::__clzll(code ^ node_code[idx - 1]);
  int prefix_next = ::__clzll(code ^ node_code[idx + 1]);

  int d = (prefix_next > prefix_prev) ? 1 : -1;
  int prefix_min = min(prefix_prev, prefix_next);

  // Find upper bound for length.

  int lmax = 128 >> 2;
  unsigned int probe;
  do
  {
    lmax <<= 2;
    probe = idx + lmax * d;
  }
  while(probe < num_leaves && ::__clzll(code ^ node_code[probe]) > prefix_min);

  // Determine length.

  int l = 0;
  for (int t = lmax >> 1; t > 0; t >>= 1)
  {
    probe = idx + (l + t) * d;
    if (probe < num_leaves && ::__clzll(code ^ node_code[probe]) > prefix_min)
        l += t;
  }
  // int j = idx + l * d;

  uint32_t jdx = idx + l * d;
  if(d < 0)
  {
      swap(idx, jdx); // make it sure that idx < jdx
  }

  return make_uint2(idx, jdx);
}

__device__ __forceinline__
uint2 determine_range(const uint64_t* node_code,
        const uint32_t num_leaves, uint32_t idx)
{
  if(idx == 0)
  {
      return make_uint2(0, num_leaves-1);
  }

  // determine direction of the range
  const uint64_t self_code = node_code[idx];
  const int L_delta = ::__clzll(self_code ^ node_code[idx-1]);
  const int R_delta = ::__clzll(self_code ^ node_code[idx+1]);
  const int d = (R_delta > L_delta) ? 1 : -1;

  // Compute upper bound for the length of the range

  // const int delta_min = min(L_delta, R_delta);
  const int delta_min = ::__clzll(self_code ^ node_code[idx-d]);
  int l_max = 2;
  int delta = -1;
  int i_tmp = idx + l_max * d;
  if(0 <= i_tmp && i_tmp < num_leaves)
  {
      delta = ::__clzll(self_code ^ node_code[i_tmp]);
  }
  while(delta > delta_min)
  {
      l_max <<= 1;
      delta = -1;
      i_tmp = idx + l_max * d;
      if(0 <= i_tmp && i_tmp < num_leaves)
      {
          delta = ::__clzll(self_code ^ node_code[i_tmp]);
      }
  }

  // Find the other end by binary search
  uint32_t l = 0;
  // int t = l_max >> 1;
  for(uint32_t t= l_max >> 1; t>0; t>>=1)
  {
    i_tmp = idx + (l + t) * d;
    delta = -1;
    if(0 <= i_tmp && i_tmp < num_leaves)
    {
      if( ::__clzll(self_code ^ node_code[i_tmp]) > delta_min )
        l += t;
    }
  }
  // while(t > 0)
  // {
  //     i_tmp = idx + (l + t) * d;
  //     delta = -1;
  //     if(0 <= i_tmp && i_tmp < num_leaves)
  //     {
  //         delta = ::__clzll(self_code ^ node_code[i_tmp]);
  //     }
  //     if(delta > delta_min)
  //     {
  //         l += t;
  //     }
  //     t >>= 1;
  // }
  uint32_t jdx = idx + l * d;
  if(d < 0)
  {
      swap(idx, jdx); // make it sure that idx < jdx
  }

#ifdef DEBUG
  printf("L_d: %d, R_d: %d, d: %d, range: [%u,%u]\n",
          L_delta, R_delta, d, idx, jdx);
#endif


  return make_uint2(idx, jdx);
}


__device__ __forceinline__
uint32_t find_split(const uint64_t* node_code, const uint32_t num_leaves,
    const uint32_t first, const uint32_t last) noexcept
{
    const uint64_t first_code = node_code[first];
    const uint64_t last_code  = node_code[last];
    if (first_code == last_code)
    {
        return (first + last) >> 1;
    }   
    const int delta_node = ::__clzll(first_code ^ last_code);

    // binary search...
    int split  = first;
    int stride = last - first;
    do
    {
        stride = (stride + 1) >> 1;
        const int middle = split + stride;
        if( middle < last && ::__clzll(first_code ^ node_code[middle]) > delta_node)
          split = middle;

    }
    while(stride > 1);

    return split;
}

__device__ __forceinline__
uint2 build_internal(const uint64_t* node_code,
        const uint32_t num_leaves, uint32_t idx)
{
  uint32_t jdx;
  if(idx != 0)
  {
    // determine direction of the range
    const uint64_t self_code = node_code[idx];
    const int L_delta = ::__clzll(self_code ^ node_code[idx-1]);
    const int R_delta = ::__clzll(self_code ^ node_code[idx+1]);
    const int d = (R_delta > L_delta) ? 1 : -1;

    // Compute upper bound for the length of the range

    const int delta_min = ::__clzll(self_code ^ node_code[idx-d]);
    int l_max = 2;
    int delta = -1;
    int i_tmp = idx + l_max * d;
    if(0 <= i_tmp && i_tmp < num_leaves)
    {
        delta = ::__clzll(self_code ^ node_code[i_tmp]);
    }
    while(delta > delta_min)
    {
        l_max <<= 1;
        delta = -1;
        i_tmp = idx + l_max * d;
        if(0 <= i_tmp && i_tmp < num_leaves)
        {
            delta = ::__clzll(self_code ^ node_code[i_tmp]);
        }
    }

    // Find the other end by binary search
    uint32_t l = 0;
    for(uint32_t t= l_max >> 1; t>0; t>>=1)
    {
      i_tmp = idx + (l + t) * d;
      delta = -1;
      if(0 <= i_tmp && i_tmp < num_leaves)
      {
        if( ::__clzll(self_code ^ node_code[i_tmp]) > delta_min )
          l += t;
      }
    }

    jdx = idx + l * d;
    if(d < 0)
    {
        swap(idx, jdx); // make it sure that idx < jdx
    }
  }
  else
  {
    jdx = num_leaves-1;
  }

  const uint64_t first_code = node_code[idx];
  const uint64_t last_code  = node_code[jdx];
  uint32_t split;
  if (first_code == last_code)
  {
    split = (idx + jdx) >> 1;
  }
  else{
    const uint32_t delta_node = ::__clzll(first_code ^ last_code);

    // binary search...
    split  = idx;
    uint32_t stride = jdx - idx;
    do
    {
        stride = (stride + 1) >> 1;
        const uint32_t middle = split + stride;
        if( middle < jdx && ::__clzll(first_code ^ node_code[middle]) > delta_node)
          split = middle;

    }
    while(stride > 1);
  }

  uint32_t lidx = split;
  uint32_t ridx = split+1;

  if(min(idx, jdx) == lidx)
  {
      lidx += num_leaves-1;
  }
  if(max(idx, jdx) == ridx)
  {
      ridx += num_leaves-1;
  }

  return make_uint2(lidx, ridx);
}


__global__
void init_nodes( node_t* n, const uint64_t* m64, const size_t num_objects )
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  // const uint32_t stride = blockDim.x * gridDim.x;

  const size_t num_internal_nodes = num_objects-1;

  if(index < num_internal_nodes) 
  {
  // for(uint32_t i = index; i<num_internal_nodes; i += stride)
  // {
    const uint2 ij  = determine_range(m64, num_objects, index);

    // find gamma
    uint32_t lidx = find_split(m64, num_objects, ij.x, ij.y);

    uint32_t ridx = lidx + 1;

    if(min(ij.x, ij.y) == lidx)
    {
        lidx += num_internal_nodes;
    }
    if(max(ij.x, ij.y) == ridx)
    {
        ridx += num_internal_nodes;
    }
    // node_listS[idx] = aux;

    n[index].left_idx = lidx;
    n[index].right_idx = ridx;
    n[lidx].parent_idx  = index;
    n[ridx].parent_idx = index;
  }

  return;

}

__global__
void init_nodes2( node_t* n, const uint64_t* node_code, const uint32_t num_leaves )
{
  const uint32_t i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i >= num_leaves-1)
    return;
    
  // const uint2 ij = build_internal(node_code, num_leaves, idx);

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

  } while( 0 <= i_tmp && i_tmp < num_leaves && delta_min < ::__clzll(self_code ^ node_code[i_tmp]) );


  // Find the other end by binary search
  uint32_t l = 0;
  for(uint32_t t= l_max >> 1; t>0; t>>=1)
  {
    i_tmp = idx + (l + t) * d;
    if( 0 <= i_tmp && i_tmp < num_leaves && delta_min < ::__clzll(self_code ^ node_code[i_tmp]))
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
  printf("ID: %u, L_d: %d, R_d: %d, d: %d, range: [%u,%u], gamma: %u, prefix: %u\n",
          i, L_delta, R_delta, d, idx, jdx, gamma, prefix_node);
#endif


  uint32_t lidx = (idx == gamma) ? (gamma + num_leaves - 1) : gamma;
  uint32_t ridx = (jdx == gamma + 1) ? (gamma + num_leaves) : (gamma + 1);

  n[i].left_idx = lidx;
  n[i].right_idx = ridx;
  n[lidx].parent_idx  = i;
  n[ridx].parent_idx = i;


  return;

}

__global__
void init_nodes3( node_t* n, const uint64_t* node_code, const uint32_t num_leaves )
{
  const uint32_t i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < num_leaves-1)
  {
    // const uint2 ij = build_internal(node_code, num_leaves, idx);

    uint32_t idx = i;
    // Choose direction.
    unsigned int code = node_code[idx];
    int prefix_prev = (i==0)? FW_S32_MIN : ::__clzll(code ^ node_code[idx - 1]);
    int prefix_next = ::__clzll(code ^ node_code[idx + 1]);

    int d = (prefix_next > prefix_prev) ? 1 : -1;
    int prefix_min = min(prefix_prev, prefix_next);

    // Find upper bound for length.

    int lmax = 128 >> 2;
    unsigned int probe;
    do
    {
      lmax <<= 2;
      probe = idx + lmax * d;
    }
    while(probe < num_leaves && ::__clzll(code ^ node_code[probe]) > prefix_min);

    // Determine length.

    int l = 0;
    for (int t = lmax >> 1; t > 0; t >>= 1)
    {
      probe = idx + (l + t) * d;
      if (probe < num_leaves && ::__clzll(code ^ node_code[probe]) > prefix_min)
          l += t;
    }
    // int j = idx + l * d;

    uint32_t jdx = idx + l * d;
    // if(d < 0)
    // {
    //     swap(idx, jdx); // make it sure that idx < jdx
    // }
  
    // const uint64_t first_code = node_code[idx];
    // const uint64_t last_code  = node_code[jdx];
    // const uint32_t delta_node = ::__clzll(first_code ^ last_code);

    // // binary search...
    // uint32_t split  = idx;
    // uint32_t stride = jdx - idx;
    // do
    // {
    //     stride = (stride + 1) >> 1;
    //     const uint32_t middle = split + stride;
    //     if( middle < jdx && delta_node < ::__clzll(first_code ^ node_code[middle]))
    //       split = middle;

    // }
    // while(stride > 1);

    int prefix_node = ::__clzll(code ^ node_code[jdx]);

    int s = 0;
    int t = l;
    do
    {
        t = (t + 1) >> 1;
        probe = idx + (s + t) * d;
        if (probe < (unsigned int)num_leaves && prefix_node < ::__clzll(code ^ node_code[probe]))
            s += t;
    }
    while (t > 1);
    int k = idx + s * d + min(d, 0);

    // Output node.

    // int lo = min(idx, jdx);
    // int hi = max(idx, jdx);
  
    uint32_t lidx = (min(idx, jdx) == k) ? (k + num_leaves - 1) : k;
    uint32_t ridx = (max(idx, jdx) == k + 1) ? (k + num_leaves) : (k + 1);

  
    n[i].left_idx = lidx;
    n[i].right_idx = ridx;
    n[lidx].parent_idx  = i;
    n[ridx].parent_idx = i;

  }

  return;

}

__global__
void init_nodes4( node_t* n, const uint64_t* node_code, const uint32_t num_leaves )
{
  const uint32_t i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i >= num_leaves-1)
    return;
    
  // const uint2 ij = build_internal(node_code, num_leaves, idx);

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

  } while( 0 <= i_tmp && i_tmp < num_leaves && delta_min < ::__clzll(self_code ^ node_code[i_tmp]) );


  // Find the other end by binary search
  uint32_t l = 0;
  for(uint32_t t= l_max >> 1; t>0; t>>=1)
  {
    i_tmp = idx + (l + t) * d;
    if( 0 <= i_tmp && i_tmp < num_leaves && delta_min < ::__clzll(self_code ^ node_code[i_tmp]))
      l += t;
  }

  uint32_t jdx = idx + l * d;
  // if(d < 0)
  // {
  //     swap(idx, jdx); // make it sure that idx < jdx
  // }

  // const uint64_t first_code = node_code[idx];
  // const uint64_t last_code  = node_code[jdx];
  const uint32_t prefix_node = ::__clzll(self_code ^ node_code[jdx]);

  // binary search...
  uint32_t split  = 0;
  uint32_t stride = l;

  do
  {
      stride = (stride + 1) >> 1;
      i_tmp = idx + (split + stride) * d;
      if( 0 <= i_tmp && i_tmp < num_leaves && prefix_node < ::__clzll(self_code ^ node_code[i_tmp]))
        split += stride;

  } while(stride > 1);
  
  const uint32_t k = idx + split * d + min(d, 0);

  uint32_t lidx = (min(idx,jdx) == k) ? (k + num_leaves - 1) : k;
  uint32_t ridx = (max(idx,jdx) == k + 1) ? (k + num_leaves) : (k + 1);

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
    while(parent != 0xFFFFFFFFu) // means idx == 0
    {
      const int old = atomicCAS(flags + parent, 0u, 1u);
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
void create_aabb_size( node_t* n, aabb_t* bb, uint32_t* flags, const point_t* p, 
      const size_t num_internal_nodes, const size_t num_nodes)
{
  uint32_t index = num_internal_nodes + threadIdx.x + blockIdx.x * blockDim.x;
  // uint32_t stride = blockDim.x * gridDim.x;
  // printf("index %u, threadIdx, %d, blockIdx %d\n", index, threadIdx.x, blockIdx.x);

  if(index < num_nodes)
  {
    uint32_t parent = n[index].parent_idx;
    while(parent != 0xFFFFFFFFu) // means idx == 0
    {
      const uint32_t old = atomicCAS(flags + parent, 0u, 1u);
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

    while(parent != 0xFFFFFFFFu && atomicCAS(flags + parent, 0, 1)) 
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

/* Saves the execution time and other markers of the algorithm  */
int save_time(std::string file_name, std::string map_name, uint32_t chunkDim, uint32_t threadsPerBlock,
  uint32_t numberOfBlocks, double build_time, int problems)
{
  // uint32_t size = pointList.size();
  std::fstream out(file_name, out.out | out.app);
  // std::ofstream out;
  // out.precision(std::numeric_limits<double>::digits10);
  // out.open(file_name);
  // Point aux;
  if(!out.is_open())
    return -1;

  out << map_name << " " << chunkDim << " " << threadsPerBlock << " ";
  out << numberOfBlocks << " " << build_time << " " << problems << std::endl;

  // out << map_name << " " << cpu_cores << " " << is_gpu_used << " ";
  // // out.precision(3);
  // out << std::defaultfloat << minRadius << " " << maxSize << " " << level << " ";
  // // out.precision(6);
  // out << std::fixed << cputree_time << " " << gputree_time << " " << time << " ";
  // // out.precision(1);
  // out << std::defaultfloat << rate << " " << check_rate << std::endl;

  out.close();
  if(out.is_open())
    return -1;
  return 0;
}



int main(int argc, char* argv[])
{
  std::string inputTXT = (argc > 1)? argv[1] : "data/INAER_2011_Alcoy.xyz";
  std::size_t N;
  aabb_t BBox;
  int problems = 0;

  node_t default_node;
  default_node.parent_idx = 0xFFFFFFFFu;
  default_node.left_idx   = 0xFFFFFFFFu;
  default_node.right_idx  = 0xFFFFFFFFu;
  default_node.numPts     = 0xFFFFFFFFu;
  default_node.min_idx    = 0xFFFFFFFFu;
  default_node.object_idx = 0xFFFFFFFFu;

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

  // #ifdef DEBUG
  //   BBox.lower.x = 0.0;
  //   BBox.lower.y = 0.0;

  //   BBox.upper.x = 10.0;
  //   BBox.upper.y = 10.0;
  // #endif

  N = (argc > 5)? static_cast<size_t>(atoi(argv[5])) : N; 

  real_v diffBox;
  diffBox.x = BBox.upper.x - BBox.lower.x;
  diffBox.y = BBox.upper.y - BBox.lower.y;
  // diffBox.z = BBox.upper.z - BBox.lower.z;

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, NUM_PROCS);

  // std::vector<point_t> ps(N);

  // if(read_pointsC(inputTXT, ps) < 0){
  //     printf("Unable to read file!\n");
  //     exit(-1);
  // }

  point_t*   point_cloud = NULL;
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
  uint32_t* flagsS;

  uint64_t* morton64S = NULL;
  node_t* node_listS = NULL;
  aabb_t* aabb_listS = NULL;

  const size_t num_objects = N;
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

  const int chunkDim = (argc > 2)? atoi(argv[2]) : 256;
  const uint32_t r_num_objects = (uint32_t)((num_objects-1)/chunkDim) + 1;
  const uint32_t r_num_internal_nodes = r_num_objects - 1;
  const uint32_t r_num_nodes = 2*r_num_objects -1;

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
  cudaMallocManaged(&flagsS, r_num_internal_nodes*sizeof(uint32_t));

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
  cudaMemPrefetchAsync(flagsS, r_num_internal_nodes*sizeof(uint32_t), deviceId);

  cudaMemPrefetchAsync(point_cloud, size_points, cudaCpuDeviceId);

#ifdef DEBUG
  for(int i=0; i<N; i++){
    point_cloud[i].x = (int)uni(mt);
    point_cloud[i].y = (int)uni(mt);
    point_cloud[i].z = (int)uni(mt);
  }
  point_cloud[9].z = -1.0;
#else
  map_file_atof_tbb_th(inputTXT.c_str(), num_objects, point_cloud);
#endif
  
  cudaMemPrefetchAsync(point_cloud, size_points, deviceId);

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


  std::cout << BBox.upper.x << "," << BBox.upper.y << " ";
  std::cout << BBox.lower.x << "," << BBox.lower.y << "\n";

  std::cout << diffBox.x << "," << diffBox.y << "\n";

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
  cudaMalloc(&d_temp_storage, temp_storage_bytes);
  // cudaMemPrefetchAsync(d_temp_storage, temp_storage_bytes, deviceId);
  
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
  // for(int i = 0; i<N; i++){
  //   std::cout << std::bitset<64>(morton64[i]) << " ";
  // }
  std::cout << "\n";
  for(int i = 0; i<N; i++){
    std::cout << indices_out[i] << " ";
  }
  std::cout << "\n";
  for(int i = 0; i<r_num_objects; i++){
    std::cout << std::bitset<64>(morton64S[i]) << " ";
  }
  std::cout << "\n";
  std::cout << "\n";
#endif

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

  // init_nodes<<<nBlocks_node_func, threadsPerBlock>>>( node_listS, morton64S, r_num_objects );
  init_nodes2<<<nBlocks_node_func, threadsPerBlock>>>( node_listS, morton64S, r_num_objects );
  // init_nodes3<<<nBlocks_node_func, threadsPerBlock>>>( node_listS, morton64S, r_num_objects );
  // init_nodes4<<<nBlocks_node_func, threadsPerBlock>>>( node_listS, morton64S, r_num_objects );

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

#ifdef DEBUG
  for(int i = 0; i<r_num_nodes; i++){
    std::cout << node_listS[i].parent_idx << ", " << node_listS[i].left_idx << ", ";
    std::cout << node_listS[i].right_idx << ", " << node_listS[i].object_idx << ", ";
    std::cout << node_listS[i].numPts << ", " << node_listS[i].min_idx << "\n";
  }
  std::cout << std::endl;
  // for(int i = 0; i<r_num_nodes; i++){
  //     std::cout << aabb_listS[i].upper.x << "," << aabb_listS[i].upper.y << " ";
  //     std::cout << aabb_listS[i].lower.x << "," << aabb_listS[i].lower.y << "\n";
  //   }
  std::cout << std::endl;
#endif
  
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
  if(asyncErr != cudaSuccess){
    printf("Error AABBs sync: %s\n", cudaGetErrorString(asyncErr));
    problems++;
  } 

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " Create AABB time elapsed: " << mytime << " ms\n";
  totaltime += mytime;
  std::cout << "  CREATION takes: " << totaltime << " ms\n";

#ifdef DEBUG
  for(int i = 0; i<r_num_nodes; i++){
    std::cout << node_listS[i].parent_idx << ", " << node_listS[i].left_idx << ", ";
    std::cout << node_listS[i].right_idx << ", " << node_listS[i].object_idx << ", ";
    std::cout << node_listS[i].numPts << ", " << node_listS[i].min_idx << "(";
    if(node_listS[i].min_idx != 0xFFFFFFFFu)
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

  if(save_time("resultados.csv", inputTXT, chunkDim, threadsPerBlock, numberOfBlocks, totaltime, problems) < 0){

    printf("Unable to create results file!\n");

  }


  cudaFree(point_cloud);
  cudaFree(morton);
  cudaFree(morton_out);
  cudaFree(indices);
  cudaFree(indices_out);
  cudaFree(d_temp_storage);

  cudaFree(flagsS);

  cudaFree(morton64S);
  cudaFree(node_listS);
  cudaFree(aabb_listS);

  


  return 0;
}
