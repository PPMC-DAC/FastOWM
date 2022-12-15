#include "lbvh.cuh"
#include <random>
#include <vector>
#include <chrono>
#include <thrust/random.h>
#include <tbb/parallel_invoke.h>
#define TBB_PREVIEW_GLOBAL_CONTROL 1

#include "tbb/global_control.h"

// #include <cub/cub.cuh>
// #include <cub/device/device_radix_sort.cuh>
#include <bitset>

using real_t = double;
using real_v = double2;
using point_t = double4;
using node_t = lbvh::detail::node;

struct aabb_t
{
    real_v upper;
    real_v lower;
};

#define TOL 1e-10
#define FW_S32_MIN  (~0x7FFFFFFF)


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
inline uint32_t morton2D(const point_t xyz /*, const real_t resolution = 1024.0f*/ ) noexcept
{
    const uint32_t xx = SeparateBy1(static_cast<uint32_t>(xyz.x));
    const uint32_t yy = SeparateBy1(static_cast<uint32_t>(xyz.y));
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

    aux.x -= BBox.lower.x;
    aux.y -= BBox.lower.y;
    aux.x /= diffBox.x;
    aux.y /= diffBox.y;

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
void create_aabb2( const node_t* n, aabb_t* bb, int* flags, 
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
      bb[parent] = merge(lbox, rbox);

      // look the next parent...
      parent = n[parent].parent_idx;
    }

  }

  return;
}

void create_aabb_cpu( const node_t* n, aabb_t* bb, const uint32_t current)
{
  if(n[current].object_idx != 0xFFFFFFFF)
  {
    return;
  }
  else
  {
    const uint32_t lidx = n[current].left_idx;
    const uint32_t ridx = n[current].right_idx;
    tbb::parallel_invoke([&]() { create_aabb_cpu( n, bb, lidx);}, [&]() { create_aabb_cpu( n, bb, ridx);} );
    // create_aabb_cpu( n, bb, lidx);
    // create_aabb_cpu( n, bb, ridx);
    const aabb_t lbox = bb[lidx];
    const aabb_t rbox = bb[ridx];
    bb[current] = merge(lbox, rbox);
  }
  return;
}


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

__device__ __inline__ int   flo     (unsigned int v)        { unsigned int r; asm("bfind.u32 %0, %1;" : "=r"(r) : "r"(v)); return r; }
__device__ __inline__ int   slct    (int a, int b, int c)   { int v; asm("slct.s32.s32 %0, %1, %2, %3;" : "=r"(v) : "r"(a), "r"(b), "r"(c)); return v; }

__device__ __inline__ int clz(unsigned int hi, unsigned int lo)
{
    int a = flo(hi);
    int b = flo(lo);
    return slct(31 - a, 63 - b, a);
}

__global__ 
void build(int4* nodeChildren, int* nodeParent, int2* nodeRange, unsigned int* nodeCounter, const unsigned int* objectCodes, int numObjects)
{
    // Pick node.

    int i = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * (blockIdx.x + gridDim.x * blockIdx.y));
    if (i >= numObjects - 1)
        return;

    // Choose direction.

    unsigned int code = objectCodes[i];
    int prefix_prev = (i == 0) ? FW_S32_MIN : clz(code ^ objectCodes[i - 1], i ^ (i - 1));
    int prefix_next = clz(code ^ objectCodes[i + 1], i ^ (i + 1));

    int d = (prefix_next > prefix_prev) ? 1 : -1;
    int prefix_min = min(prefix_prev, prefix_next);

    // Find upper bound for length.

    int lmax = 128 >> 2;
    unsigned int probe;
    do
    {
        lmax <<= 2;
        probe = i + lmax * d;
    }
    while (probe < (unsigned int)numObjects && clz(code ^ objectCodes[probe], i ^ probe) > prefix_min);

    // Determine length.

    int l = 0;
    for (int t = lmax >> 1; t > 0; t >>= 1)
    {
        probe = i + (l + t) * d;
        if (probe < (unsigned int)numObjects & clz(code ^ objectCodes[probe], i ^ probe) > prefix_min)
            l += t;
    }
    int j = i + l * d;
    int prefix_node = clz(code ^ objectCodes[j], i ^ j);

    // Find split point.

    int s = 0;
    int t = l;
    do
    {
        t = (t + 1) >> 1;
        probe = i + (s + t) * d;
        if (probe < (unsigned int)numObjects && clz(code ^ objectCodes[probe], i ^ probe) > prefix_node)
            s += t;
    }
    while (t > 1);
    int k = i + s * d + min(d, 0);

    // Output node.

    int lo = min(i, j);
    int hi = max(i, j);

    int4 c;
    c.x = (lo == k) ? (k + numObjects - 1) : k;
    c.y = (hi == k + 1) ? (k + numObjects) : (k + 1);
    c.z = k;
    c.w = hi;

    nodeChildren[i] = c;
    nodeParent[c.x] = i;
    nodeParent[c.y] = i;
    nodeRange[i] = make_int2(lo, hi + 1);
    nodeCounter[i] = 0;

    if (i == 0)
        nodeParent[i] = 0;
}



int main(int argc, char* argv[])
{
  std::string inputTXT = (argc > 1)? argv[1] : "INAER_2011_Alcoy.xyz";
  std::size_t N;
  aabb_t BBox;

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

  N = (argc > 2)? static_cast<size_t>(atoi(argv[2])) : N; 

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, 6);

  real_v diffBox;
  diffBox.x = BBox.upper.x - BBox.lower.x;
  diffBox.y = BBox.upper.y - BBox.lower.y;
  // diffBox.z = BBox.upper.z - BBox.lower.z;


  std::vector<point_t> ps(N);

  std::cout << N << " points\n";

  if(read_pointsC(inputTXT, ps) < 0){
      printf("Unable to read file!\n");
      exit(-1);
  }

  point_t*   point_cloud = NULL;
  uint32_t*  morton = NULL;
  uint64_t*  morton64 = NULL;
  uint32_t*  morton_out = NULL;
  uint32_t*  indices = NULL;
  uint32_t*  indices_out = NULL;
  void*      d_temp_storage = NULL;

  node_t* node_list = NULL;
  node_t* aux_node_list = NULL;
  aabb_t* aabb_list = NULL;
  aabb_t* aux_aabb_list = NULL;

  int* flags;

  const size_t num_objects = N;
  const size_t num_internal_nodes = num_objects - 1;
  const size_t num_nodes = 2*num_objects - 1; /*Numero de Nodos*/ 

  int deviceId;
  int numberOfSMs;

  const size_t size_points = num_objects * sizeof(point_t);
  const size_t size_morton = num_objects * sizeof(uint32_t);
  const size_t size_morton64 = num_objects * sizeof(uint64_t);
  size_t temp_storage_bytes = 0;
  
  const size_t size_nodes = num_nodes * sizeof(node_t);
  const size_t size_aabbs = num_nodes * sizeof(aabb_t);

  cudaError_t mortonError;
  cudaError_t asyncErr;

  cudaGetDevice(&deviceId);
  cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);
  printf("Device ID: %d\tNumber of SMs: %d\n", deviceId, numberOfSMs);

  size_t threadsPerBlock = 256;
  size_t numberOfBlocks = 32 * numberOfSMs;

  cudaMallocManaged(&point_cloud, size_points);
  cudaMallocManaged(&morton, size_morton);
  cudaMallocManaged(&morton64, size_morton64);
  cudaMallocManaged(&morton_out, size_morton);
  cudaMallocManaged(&indices, size_morton);
  cudaMallocManaged(&indices_out, size_morton);

  cudaMallocManaged(&node_list, size_nodes);
  aux_node_list = &node_list[num_internal_nodes];
  cudaMallocManaged(&aabb_list, size_aabbs);
  aux_aabb_list = &aabb_list[num_internal_nodes];

  cudaMallocManaged(&flags, N*sizeof(int));

  cudaMemPrefetchAsync(node_list, size_nodes, deviceId);
  cudaMemPrefetchAsync(aabb_list, size_aabbs, deviceId);
  cudaMemPrefetchAsync(morton, size_morton, deviceId);
  cudaMemPrefetchAsync(indices, size_morton, deviceId);
  cudaMemPrefetchAsync(morton_out, size_morton, deviceId);
  cudaMemPrefetchAsync(indices_out, size_morton, deviceId);
  cudaMemPrefetchAsync(morton64, size_morton64, deviceId);
  cudaMemPrefetchAsync(flags, num_objects*sizeof(int), deviceId);

  cudaMemPrefetchAsync(point_cloud, size_points, cudaCpuDeviceId);

  int i = 0;
  for(auto& p : ps){
    point_cloud[i].x = p.x;
    point_cloud[i].y = p.y;
    point_cloud[i].z = p.z;
    i++;
  }

  cudaMemPrefetchAsync(point_cloud, size_points, deviceId);

  // int i = 0;
  // for(auto& p : ps){
  //   point_cloud[i].x = (int)uni(mt);
  //   point_cloud[i].y = (int)uni(mt);
  //   point_cloud[i].z = (int)uni(mt);
  //   i++;
  // }

  node_t default_node;
  default_node.parent_idx = 0xFFFFFFFF;
  default_node.left_idx   = 0xFFFFFFFF;
  default_node.right_idx  = 0xFFFFFFFF;
  default_node.object_idx = 0xFFFFFFFF;

  const auto inf = std::numeric_limits<double>::infinity();
  aabb_t default_aabb;
  default_aabb.upper.x = -inf; default_aabb.lower.x = inf;
  default_aabb.upper.y = -inf; default_aabb.lower.y = inf;
  // default_aabb.upper.z = -inf; default_aabb.lower.z = inf;

  // cudaMemPrefetchAsync(node_list, size_nodes, deviceId);
  // cudaMemPrefetchAsync(aabb_list, size_aabbs, deviceId);

  init_struct<<<numberOfBlocks, threadsPerBlock>>>( node_list, aabb_list, num_nodes, default_node, default_aabb);
  // init_struct<<<numberOfBlocks, threadsPerBlock>>>( aabb_list, num_internal_nodes, default_aabb);
  // init_struct<<<numberOfBlocks, threadsPerBlock>>>( node_list, num_internal_nodes, default_node);

  mortonError = cudaGetLastError();
  if(mortonError != cudaSuccess) printf("Error INIT NODES y ABBBs: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error INIT NODES y ABBBs sync: %s\n", cudaGetErrorString(asyncErr));


  // cudaMemPrefetchAsync(point_cloud, size_points, deviceId);
  // cudaMemPrefetchAsync(morton, size_morton, deviceId);
  // cudaMemPrefetchAsync(indices, size_morton, deviceId);

  double totaltime = 0;
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();

  get_morton_code<<<numberOfBlocks, threadsPerBlock>>>( point_cloud, BBox, diffBox, morton, indices, num_objects);

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error MORTON: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error MORTON sync: %s\n", cudaGetErrorString(asyncErr));

  // double mytime = cast_t(tempo_t::now() - i_start).count();
  // std::cout << "  MORTON time elapsed: " << mytime << " ms\n";
  // double totaltime = 0;

#ifdef DEBUG
  for(int i = 0; i<N; i++){
    std::cout << point_cloud[i].x << ",";
    std::cout << point_cloud[i].y << " ";
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


  // i_start = tempo_t::now();

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

  // mytime = cast_t(tempo_t::now() - i_start).count();
  // std::cout << "  SORT time elapsed: " << mytime << " ms\n";
  // totaltime += mytime;

#ifdef DEBUG
  for(int i = 0; i<N; i++){
    std::cout << std::bitset<32>(morton_out[i]) << " ";
  }
  std::cout << "\n";
  for(int i = 0; i<N; i++){
    std::cout << indices_out[i] << " ";
  }
  std::cout << "\n";
#endif

  // cudaMemPrefetchAsync(morton64, size_morton64, deviceId);

  // i_start = tempo_t::now();

  check_morton<<<numberOfBlocks, threadsPerBlock>>>( morton_out, indices_out, morton64, num_objects );

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error CHECK: %s\n", cudaGetErrorString(mortonError));

  // asyncErr = cudaDeviceSynchronize();
  // if(asyncErr != cudaSuccess) printf("Error CHECK sync: %s\n", cudaGetErrorString(asyncErr));

  // mytime = cast_t(tempo_t::now() - i_start).count();
  // std::cout << "  CHECK time elapsed: " << mytime << " ms\n";
  // totaltime += mytime;

#ifdef DEBUG
  for(int i = 0; i<N; i++){
    std::cout << std::bitset<64>(morton64[i]) << " ";
  }
  std::cout << "\n";
  for(int i = 0; i<N; i++){
    std::cout << indices_out[i] << " ";
  }
  std::cout << "\n";
#endif


  // cudaMemPrefetchAsync(aux_node_list, num_objects*sizeof(node_t), deviceId);
  // cudaMemPrefetchAsync(aux_aabb_list, num_objects*sizeof(aabb_t), deviceId);

  // i_start = tempo_t::now();

  // cudaMemsetAsync( node_list, 0xFFFFFFFF, size_nodes*0.5 , (cudaStream_t)1);

  init_leafs<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, aux_aabb_list, point_cloud, indices_out, num_objects );

  // init_leafs<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, indices, N );

  // init_leafs2<<<numberOfBlocks, threadsPerBlock>>>( node_list, aabb_list, point_cloud, indices, num_internal_nodes, num_nodes );

  // init_indices<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, indices, N );

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error LEAFS: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error LEAFS sync: %s\n", cudaGetErrorString(asyncErr));

  // mytime = cast_t(tempo_t::now() - i_start).count();
  // std::cout << "  LEAFS and AABBs time elapsed: " << mytime << " ms\n";
  // totaltime += mytime;

#ifdef DEBUG
  for(int i = 0; i<num_nodes; i++){
      std::cout << node_list[i].parent_idx << ", " << node_list[i].left_idx << ", ";
      std::cout << node_list[i].right_idx << ", " << node_list[i].object_idx << "\n";
  }
  std::cout << std::endl;

  for(int i = 0; i<num_nodes; i++){
      std::cout << aabb_list[i].upper.x << "," << aabb_list[i].upper.y << " ";
      std::cout << aabb_list[i].lower.x << "," << aabb_list[i].lower.y << "\n";
  }
  std::cout << std::endl;
#endif

  // cudaMemPrefetchAsync(node_list, size_nodes, deviceId);
  // cudaMemPrefetchAsync(morton64, size_morton64, deviceId);

  // i_start = tempo_t::now();

  // init_nodes<<<numberOfBlocks, threadsPerBlock>>>( node_list, morton64, N );

  thrust::for_each(thrust::device,
      thrust::make_counting_iterator<uint32_t>(0),
      thrust::make_counting_iterator<uint32_t>(num_internal_nodes),
      [node_list, morton64, num_objects] __device__ (const uint32_t idx)
      {
        // node_t aux;

        const uint2 ij  = lbvh::detail::determine_range(morton64, num_objects, idx);
        // const uint32_t gamma = lbvh::detail::find_split(morton64, num_objects, ij.x, ij.y);
        uint32_t left_idx = lbvh::detail::find_split(morton64, num_objects, ij.x, ij.y);

#ifdef DEBUG
        printf("%u, %u left_idx: %d\n",ij.x, ij.y, left_idx);
#endif

        // aux.object_idx = 0xFFFFFFFF;

        // node_list[idx].left_idx  = gamma;
        // node_list[idx].right_idx = gamma + 1;
        uint32_t right_idx = left_idx + 1;

        if(thrust::min(ij.x, ij.y) == left_idx)
        {
            left_idx += num_objects - 1;
        }
        if(thrust::max(ij.x, ij.y) == right_idx)
        {
            right_idx += num_objects - 1;
        }
        // node_list[idx] = aux;

        node_list[idx].left_idx = left_idx;
        node_list[idx].right_idx = right_idx;
        node_list[left_idx].parent_idx  = idx;
        node_list[right_idx].parent_idx = idx;
        return;
      });

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error NODES: %s\n", cudaGetErrorString(mortonError));

  // asyncErr = cudaDeviceSynchronize();
  // if(asyncErr != cudaSuccess) printf("Error NODES sync: %s\n", cudaGetErrorString(asyncErr));

  // mytime = cast_t(tempo_t::now() - i_start).count();
  // std::cout << "  INTERNAL NODES time elapsed: " << mytime << " ms\n";
  // totaltime += mytime;

#ifdef DEBUG
  for(int i = 0; i<num_nodes; i++){
      std::cout << node_list[i].parent_idx << ", " << node_list[i].left_idx << ", ";
      std::cout << node_list[i].right_idx << ", " << node_list[i].object_idx << "\n";
  }
  std::cout << std::endl;
#endif

  // cudaMemPrefetchAsync(flags, num_objects*sizeof(int), deviceId);
  // cudaMemPrefetchAsync(aabb_list, size_aabbs, deviceId);

  // i_start = tempo_t::now();

  uint32_t nBlocks_aabb_func = (uint32_t)(N/threadsPerBlock) + 1;

  create_aabb<<<nBlocks_aabb_func, threadsPerBlock>>>( node_list, aabb_list, flags, num_internal_nodes, num_nodes );
  // create_aabb2<<<numberOfBlocks, threadsPerBlock>>>( node_list, aabb_list, flags, num_internal_nodes, num_nodes );
  // create_aabb_cpu(node_list, aabb_list, 0);

  // thrust::for_each(thrust::device,
  //     thrust::make_counting_iterator<uint32_t>(num_internal_nodes),
  //     thrust::make_counting_iterator<uint32_t>(num_nodes),
  //     [node_list, aabb_list, flags] __device__ (const uint32_t idx)
  //     {
  //         uint32_t parent = node_list[idx].parent_idx;
  //         while(parent != 0xFFFFFFFF) // means idx == 0
  //         {
  //             const int old = atomicCAS(flags + parent, 0, 1);
  //             if(old == 0)
  //             {
  //                 // this is the first thread entered here.
  //                 // wait the other thread from the other child node.
  //                 return;
  //             }
  //             assert(old == 1);
  //             // here, the flag has already been 1. it means that this
  //             // thread is the 2nd thread. merge AABB of both childlen.

  //             const uint32_t lidx = node_list[parent].left_idx;
  //             const uint32_t ridx = node_list[parent].right_idx;
  //             const aabb_t lbox = aabb_list[lidx];
  //             const aabb_t rbox = aabb_list[ridx];
  //             aabb_list[parent] = lbvh::merge(lbox, rbox);

  //             // look the next parent...
  //             parent = node_list[parent].parent_idx;
  //         }
  //         return;
  //     });

  // mortonError = cudaGetLastError();
  // if(mortonError != cudaSuccess) printf("Error AABBs: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error AABBs sync: %s\n", cudaGetErrorString(asyncErr));

  double mytime = cast_t(tempo_t::now() - i_start).count();
  // std::cout << "  Create AABB time elapsed: " << mytime << " ms\n";
  totaltime += mytime;
  std::cout << "  CREATION takes: " << totaltime << " ms\n";

#ifdef DEBUG
  for(int i = 0; i<num_nodes; i++){
      std::cout << aabb_list[i].upper.x << "," << aabb_list[i].upper.y << " ";
      std::cout << aabb_list[i].lower.x << "," << aabb_list[i].lower.y << "\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<N; i++){
    std::cout << morton64[i] << " ";
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

  cudaFree(point_cloud);
  cudaFree(morton);
  cudaFree(morton64);
  cudaFree(morton_out);
  cudaFree(indices);
  cudaFree(indices_out);
  cudaFree(d_temp_storage);

  cudaFree(node_list);
  aux_node_list = NULL;
  cudaFree(aabb_list);
  aux_aabb_list = NULL;

  cudaFree(flags);

  // mytime = cast_t(tempo_t::now() - i_start).count();
  // double cpu_tree_time = mytime/1e3;

  // std::cout << "  CREATION time elapsed: " << cpu_tree_time << " s\n";

  // const auto bvh_dev = bvh.get_device_repr();

  // std::cout << "testing query_device:overlap ...\n";
  // thrust::for_each(thrust::device,
  //     thrust::make_counting_iterator<std::size_t>(0),
  //     thrust::make_counting_iterator<std::size_t>(N),
  //     [bvh_dev] __device__ (std::size_t idx) {
  //         unsigned int buffer[10];
  //         const auto self = bvh_dev.objects[idx];
  //         const real_t  dr = 0.1f;
  //         for(std::size_t i=1; i<10; ++i)
  //         {
  //             for(unsigned int j=0; j<10; ++j)
  //             {
  //                 buffer[j] = 0xFFFFFFFF;
  //             }
  //             const real_t r = dr * i;
  //             lbvh::aabb<real_t> query_box;
  //             query_box.lower = make_real_v(self.x-r, self.y-r, self.z-r, 0);
  //             query_box.upper = make_real_v(self.x+r, self.y+r, self.z+r, 0);
  //             const auto num_found = lbvh::query_device(
  //                     bvh_dev, lbvh::overlaps(query_box), buffer, 10);

  //             for(unsigned int j=0; j<10; ++j)
  //             {
  //                 const auto jdx    = buffer[j];
  //                 if(j >= num_found)
  //                 {
  //                     assert(jdx == 0xFFFFFFFF);
  //                     continue;
  //                 }
  //                 else
  //                 {
  //                     assert(jdx != 0xFFFFFFFF);
  //                     assert(jdx < bvh_dev.num_objects);
  //                 }
  //                 const auto other  = bvh_dev.objects[jdx];
  //                 assert(fabsf(self.x - other.x) < r); // check coordinates
  //                 assert(fabsf(self.y - other.y) < r); // are in the range
  //                 assert(fabsf(self.z - other.z) < r); // of query box
  //             }
  //         }
  //         return ;
  //     });

  // std::cout << "testing query_device:nearest_neighbor ...\n";
  // thrust::for_each(thrust::device,
  //     thrust::make_counting_iterator<unsigned int>(0),
  //     thrust::make_counting_iterator<unsigned int>(N),
  //     [bvh_dev] __device__ (const unsigned int idx) {
  //         const auto self = bvh_dev.objects[idx];
  //         const auto nest = lbvh::query_device(bvh_dev, lbvh::nearest(self),
  //                                              distance_calculator());
  //         assert(nest.first != 0xFFFFFFFF);
  //         const auto other   = bvh_dev.objects[nest.first];
  //         // of course, the nearest object is itself.
  //         assert(nest.second == 0.0f);
  //         assert(self.x == other.x);
  //         assert(self.y == other.y);
  //         assert(self.z == other.z);
  //         return ;
  //    });

  // thrust::device_vector<real_v> random_points(N);
  // thrust::transform(
  //     thrust::make_counting_iterator<unsigned int>(0),
  //     thrust::make_counting_iterator<unsigned int>(N),
  //     random_points.begin(), [] __device__(const unsigned int idx) {
  //         thrust::default_random_engine rand;
  //         thrust::uniform_real_distribution<real_t> uni(0.0f, 1.0f);
  //         rand.discard(idx);
  //         const real_t x = uni(rand);
  //         const real_t y = uni(rand);
  //         const real_t z = uni(rand);
  //         return make_real_v(x, y, z, 0);
  //     });

  // thrust::for_each(random_points.begin(), random_points.end(),
  //     [bvh_dev] __device__ (const real_v pos) {
  //         const auto calc = distance_calculator();
  //         const auto nest = lbvh::query_device(bvh_dev, lbvh::nearest(pos), calc);
  //         assert(nest.first != 0xFFFFFFFF);

  //         for(unsigned int i=0; i<bvh_dev.num_objects; ++i)
  //         {
  //             const auto dist = calc(bvh_dev.objects[i], pos);
  //             if(i == nest.first)
  //             {
  //                 assert(dist == nest.second);
  //             }
  //             else
  //             {
  //                 assert(dist >= nest.second);
  //             }
  //         }
  //         return ;
  //     });

  // std::cout << "testing query_host:overlap ...\n";
  // {
  //     for(std::size_t i=0; i<10; ++i)
  //     {
  //         const auto self = bvh.objects_host()[i];
  //         const real_t dr = 0.1f;
  //         for(unsigned int cnt=1; cnt<10; ++cnt)
  //         {
  //             const real_t r = dr * cnt;
  //             lbvh::aabb<real_t> query_box;
  //             query_box.lower = make_real_v(self.x-r, self.y-r, self.z-r, 0);
  //             query_box.upper = make_real_v(self.x+r, self.y+r, self.z+r, 0);

  //             std::vector<std::size_t> buffer;
  //             const auto num_found = lbvh::query_host(bvh,
  //                     lbvh::overlaps(query_box), std::back_inserter(buffer));

  //             for(unsigned int jdx : buffer)
  //             {
  //                 assert(jdx < bvh.objects_host().size());

  //                 const auto other  = bvh.objects_host()[jdx];
  //                 assert(fabsf(self.x - other.x) < r); // check coordinates
  //                 assert(fabsf(self.y - other.y) < r); // are in the range
  //                 assert(fabsf(self.z - other.z) < r); // of query box
  //             }
  //         }
  //     }
  // }

  return 0;
}
