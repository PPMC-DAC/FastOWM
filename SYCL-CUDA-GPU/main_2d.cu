#include "lbvh.cuh"
#include <random>
#include <vector>
#include <chrono>
#include <thrust/random.h>
// #include <cub/cub.cuh>
// #include <cub/device/device_radix_sort.cuh>
#include <bitset>

using real_t = double;
using real_s = double2;
using points_t = double4;

struct aabb_t
{
    real_s upper;
    real_s lower;
};

using node_t = lbvh::detail::node;
// using aabb_t = lbvh::aabb<real_t>;

struct aabb_getter
{
  /* aqui destaca que el limite superior e inferior que le da a las cajas es el mismo,
  es decir, la caja es un solo punto */
  __device__
  aabb_t operator()(const real_s f) const noexcept
  {
      aabb_t retval;
      retval.upper = f;
      retval.lower = f;
      return retval;
  }
};
struct distance_calculator
{
    __device__
    real_t operator()(const real_s point, const real_s object) const noexcept
    {
        return (point.x - object.x) * (point.x - object.x) +
               (point.y - object.y) * (point.y - object.y);
    }
};

int read_pointsC(std::string file_name, std::vector<points_t>& point_cloud)
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
inline uint32_t morton2D(points_t xyz /*, const real_t resolution = 1024.0f*/ ) noexcept
{
    const uint32_t xx = SeparateBy1(static_cast<uint32_t>(xyz.x));
    const uint32_t yy = SeparateBy1(static_cast<uint32_t>(xyz.y));
    return xx * 2 + yy;
}

// // Expands a 10-bit integer into 30 bits
// // by inserting 2 zeros after each bit.
// __device__ __host__
// inline std::uint32_t expandBits(std::uint32_t v) noexcept
// {
//     v = (v * 0x00010001u) & 0xFF0000FFu;
//     v = (v * 0x00000101u) & 0x0F00F00Fu;
//     v = (v * 0x00000011u) & 0xC30C30C3u;
//     v = (v * 0x00000005u) & 0x49249249u;
//     return v;
// }

// // Calculates a 30-bit Morton code for the
// // given 3D point located within the unit cube [0,1].
// __device__
// inline std::uint32_t morton3D(real_s xyz, real_t resolution = 1024.0f) noexcept
// {
//     xyz.x = ::fminf(::fmaxf(xyz.x * resolution, 0.0f), resolution - 1.0f);
//     xyz.y = ::fminf(::fmaxf(xyz.y * resolution, 0.0f), resolution - 1.0f);
//     xyz.z = ::fminf(::fmaxf(xyz.z * resolution, 0.0f), resolution - 1.0f);
//     const std::uint32_t xx = expandBits(static_cast<std::uint32_t>(xyz.x));
//     const std::uint32_t yy = expandBits(static_cast<std::uint32_t>(xyz.y));
//     const std::uint32_t zz = expandBits(static_cast<std::uint32_t>(xyz.z));
//     return xx * 4 + yy * 2 + zz;
// }

__global__ 
void get_morton_code(points_t* points, aabb_t BBox, real_s diffBox, 
                        uint32_t* morton, uint32_t* indices, const size_t N)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for(int i = index; i<N; i += stride){

    indices[i] = i;

    points_t aux = points[i];

    aux.x -= BBox.lower.x;
    aux.y -= BBox.lower.y;
    aux.x /= diffBox.x;
    aux.y /= diffBox.y;

    morton[i] = morton2D(aux);

  }

  return;
}

// __global__
// void check_morton(uint32_t* m, uint32_t* idx, uint64_t* om64, const size_t N)
// {
//   int index = threadIdx.x + blockIdx.x * blockDim.x;
//   int stride = blockDim.x * gridDim.x;
//   uint64_t m64;

//   for(int i = index; i<N; i += stride)
//   {
//     m64 = m[i];
//     m64 <<= 32;
//     m64 |= idx[i];
//     om64[i] = m64;
//   }

//   return;
// }

// __global__
// void init_leafs(node_t* n, aabb_t* bb, real_s* p, uint32_t* idx, const size_t N)
// {
//   int index = threadIdx.x + blockIdx.x * blockDim.x;
//   int stride = blockDim.x * gridDim.x;

//   node_t auxn;
//   auxn.parent_idx = 0xFFFFFFFF;
//   auxn.left_idx   = 0xFFFFFFFF;
//   auxn.right_idx  = 0xFFFFFFFF;
//   // aux.object_idx = 0xFFFFFFFF;

//   for(int i = index; i<N; i += stride)
//   {
//     uint32_t id = idx[i]; // obtengo el índice del punto en esa posición

//     real_s auxp = p[id]; // obtengo el punto

//     aabb_t auxbb; // creo un BB que es solo el punto
//     auxbb.upper = auxp;
//     auxbb.lower = auxp;

//     auxn.object_idx = id; // seteo el índice del nodo

//     bb[i] = auxbb;
//     n[i] = auxn;
//   }

//   return;
// }

// __global__
// void init_leafs2(node_t* n, uint32_t* idx, const size_t N, const size_t NN)
// {
//   int index = threadIdx.x + blockIdx.x * blockDim.x;
//   int stride = blockDim.x * gridDim.x;

//   node_t aux;
//   aux.parent_idx = 0xFFFFFFFF;
//   aux.left_idx   = 0xFFFFFFFF;
//   aux.right_idx  = 0xFFFFFFFF;
//   aux.object_idx = 0xFFFFFFFF;


//   for(int i = index; i<NN; i += stride)
//   {
//     n[i] = aux;
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

// __device__
// inline int common_upper_bits(const uint64_t lhs, const uint64_t rhs) noexcept
// {
//     return ::__clzll(lhs ^ rhs);
// }

// template<typename UInt>
// __device__
// inline uint2 determine_range(UInt const* node_code,
//         const uint32_t num_leaves, uint32_t idx)
// {
//     if(idx == 0)
//     {
//         return make_uint2(0, num_leaves-1);
//     }

//     // determine direction of the range
//     const UInt self_code = node_code[idx];
//     const int L_delta = common_upper_bits(self_code, node_code[idx-1]);
//     const int R_delta = common_upper_bits(self_code, node_code[idx+1]);
//     const int d = (R_delta > L_delta) ? 1 : -1;

//     // Compute upper bound for the length of the range

//     const int delta_min = min(L_delta, R_delta);
//     int l_max = 2;
//     int delta = -1;
//     int i_tmp = idx + d * l_max;
//     if(0 <= i_tmp && i_tmp < num_leaves)
//     {
//         delta = common_upper_bits(self_code, node_code[i_tmp]);
//     }
//     while(delta > delta_min)
//     {
//         l_max <<= 1;
//         i_tmp = idx + d * l_max;
//         delta = -1;
//         if(0 <= i_tmp && i_tmp < num_leaves)
//         {
//             delta = common_upper_bits(self_code, node_code[i_tmp]);
//         }
//     }

//     // Find the other end by binary search
//     int l = 0;
//     int t = l_max >> 1;
//     while(t > 0)
//     {
//         i_tmp = idx + (l + t) * d;
//         delta = -1;
//         if(0 <= i_tmp && i_tmp < num_leaves)
//         {
//             delta = common_upper_bits(self_code, node_code[i_tmp]);
//         }
//         if(delta > delta_min)
//         {
//             l += t;
//         }
//         t >>= 1;
//     }
//     uint32_t jdx = idx + l * d;
//     if(d < 0)
//     {
//         // thrust::swap(idx, jdx); // make it sure that idx < jdx
//         uint32_t aux = idx;
//         idx = jdx;
//         jdx = aux;
//     }
//     return make_uint2(idx, jdx);
// }

// __device__
// inline uint32_t find_split(const uint64_t* m64,
//     const uint32_t first, const uint32_t last) noexcept
// {
//   // Identical Morton codes => split the range in the middle.

//   const uint64_t first_code = m64[first];
//   const uint64_t last_code  = m64[last];

//   if (first_code == last_code)
//     return (first + last) >> 1;

//   // Calculate the number of highest bits that are the same
//   // for all objects, using the count-leading-zeros intrinsic.    
  
//   // const int delta_node = common_upper_bits(first_code, last_code);
//   const int commonPrefix = __clzll(first_code ^ last_code);

//   // Use binary search to find where the next bit differs.
//   // Specifically, we are looking for the highest object that
//   // shares more than commonPrefix bits with the first one.

//   int split  = first; // initial guess
//   int step = last - first;

//   do
//   {
//     step = (step + 1) >> 1; // exponential decrease
//     const int newSplit = split + step; // proposed new position

//     if (newSplit < last)
//     {
//       const uint64_t splitCode = m64[newSplit];
//       const int splitPrefix =  __clzll(first_code ^ splitCode);

//       if (splitPrefix > commonPrefix)
//       {
//           split = newSplit; // accept proposal
//       }
//     }
//   }
//   while(step > 1);

//   return split;
// }

// __global__
// void init_nodes( node_t* n, const uint64_t* m64, const size_t N )
// {
//   uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
//   uint32_t stride = blockDim.x * gridDim.x;

//   node_t current;
//   current.parent_idx = 0xFFFFFFFF;
//   current.left_idx   = 0xFFFFFFFF;
//   current.right_idx  = 0xFFFFFFFF;
//   current.object_idx = 0xFFFFFFFF;

//   for(uint32_t i = index; i<N; i += stride)
//   {
//     const uint2 ij  = determine_range(m64, N, i);
//     const int gamma = find_split(m64, ij.x, ij.y);

//     current.left_idx  = gamma;
//     current.right_idx = gamma + 1;
//     if(min(ij.x, ij.y) == gamma)
//     {
//         current.left_idx += N - 1;
//     }
//     if(max(ij.x, ij.y) == gamma + 1)
//     {
//         current.right_idx += N - 1;
//     }
//     n[current.left_idx].parent_idx  = i;
//     n[current.right_idx].parent_idx = i;

//     n[i] = current;
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


// __global__
// void create_aabb(node_t* n, aabb_t* bb, int* flags, const size_t N, const size_t NN)
// {
//   uint32_t index = N + threadIdx.x + blockIdx.x * blockDim.x;
//   uint32_t stride = blockDim.x * gridDim.x;

//   for(uint32_t i = index; i<NN; i += stride)
//   {
//     uint32_t parent = n[i].parent_idx;
//     while(parent != 0xFFFFFFFF) // means idx == 0
//     {
//         const int old = atomicCAS(flags + parent, 0, 1);
//         if(old == 0)
//         {
//             // this is the first thread entered here.
//             // wait the other thread from the other child node.
//             return;
//         }
//         assert(old == 1);
//         // here, the flag has already been 1. it means that this
//         // thread is the 2nd thread. merge AABB of both childlen.

//         const uint32_t lidx = n[parent].left_idx;
//         const uint32_t ridx = n[parent].right_idx;
//         const aabb_t lbox = bb[lidx];
//         const aabb_t rbox = bb[ridx];
//         bb[parent] = lbvh::merge(lbox, rbox);

//         // look the next parent...
//         parent = n[parent].parent_idx;
//     }
//     return;
//   }

//   return;
// }


int main(int argc, char* argv[])
{
  std::string inputTXT = (argc > 1)? argv[1] : "INAER_2011_Alcoy.xyz";
  std::size_t N;
  aabb_t BBox;

  // const std::size_t N = (argc > 1)? atoi(argv[1]) : 10;

  std::mt19937 mt(123456789);
  // std::uniform_real_distribution<real_t> uni(0.0, 1.0);
  std::uniform_real_distribution<real_t> uni(3000.0, 5000.0);

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

  real_s diffBox;
  diffBox.x = BBox.upper.x - BBox.lower.x;
  diffBox.y = BBox.upper.y - BBox.lower.y;
  // diffBox.z = BBox.upper.z - BBox.lower.z;


  std::vector<points_t> ps(N);

  std::cout << N << " points\n";

  if(read_pointsC(inputTXT, ps) < 0){
      printf("Unable to read file!\n");
      exit(-1);
  }

  // for(auto& p : ps)
  // {
  //     p.x = uni(mt);
  //     p.y = uni(mt);
  //     p.z = uni(mt);
  // }


  points_t*  point_cloud = NULL;
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

  int deviceId;
  int numberOfSMs;
  const size_t size = N * sizeof(points_t);
  const size_t size_morton = N * sizeof(uint32_t);
  const size_t size_morton64 = N * sizeof(uint64_t);
  size_t temp_storage_bytes = 0;

  const size_t NN = 2*N; /*Numero de Nodos*/ 
  const size_t size_nodes = NN * sizeof(node_t);
  const size_t size_aabbs = NN * sizeof(aabb_t);

  cudaError_t mortonError;
  cudaError_t asyncErr;

  cudaGetDevice(&deviceId);
  cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);
  printf("Device ID: %d\tNumber of SMs: %d\n", deviceId, numberOfSMs);

  size_t threadsPerBlock = 256;
  size_t numberOfBlocks = 32 * numberOfSMs;

  cudaMallocManaged(&point_cloud, size);
  cudaMallocManaged(&morton, size_morton);
  cudaMallocManaged(&morton64, size_morton64);
  cudaMallocManaged(&morton_out, size_morton);
  cudaMallocManaged(&indices, size_morton);
  cudaMallocManaged(&indices_out, size_morton);

  cudaMallocManaged(&node_list, size_nodes);
  aux_node_list = &node_list[N];
  cudaMallocManaged(&aabb_list, size_aabbs);
  aux_aabb_list = &aabb_list[N];

  cudaMallocManaged(&flags, N*sizeof(int));

  cudaMemPrefetchAsync(point_cloud, size, cudaCpuDeviceId);

  // int i = 0;
  // for(auto& p : ps){
  //   point_cloud[i].x = i;
  //   point_cloud[i].y = i;
  //   point_cloud[i].z = i;
  //   i++;
  // }

  int i = 0;
  for(auto& p : ps){
    point_cloud[i].x = p.x;
    point_cloud[i].y = p.y;
    point_cloud[i].z = p.z;
    i++;
  }

  cudaMemPrefetchAsync(point_cloud, size, deviceId);
  cudaMemPrefetchAsync(morton, size_morton, deviceId);
  cudaMemPrefetchAsync(indices, size_morton, deviceId);

  std::chrono::time_point<tempo_t> i_start = tempo_t::now();

  get_morton_code<<<numberOfBlocks, threadsPerBlock>>>( point_cloud, BBox, diffBox, morton, indices, N);

  mortonError = cudaGetLastError();
  if(mortonError != cudaSuccess) printf("Error MORTON: %s\n", cudaGetErrorString(mortonError));

  asyncErr = cudaDeviceSynchronize();
  if(asyncErr != cudaSuccess) printf("Error MORTON sync: %s\n", cudaGetErrorString(asyncErr));

  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << "  MORTON time elapsed: " << mytime << " ms\n";
  double totaltime = mytime;

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

//   // cudaMemPrefetchAsync(morton, size_morton, deviceId);
//   cudaMemPrefetchAsync(morton_out, size_morton, deviceId);
//   // cudaMemPrefetchAsync(indices, size_morton, deviceId);
//   cudaMemPrefetchAsync(indices_out, size_morton, deviceId);


//   i_start = tempo_t::now();

//   /* Determine temporary device storage requirements; segun las especificaciones
//   si el puntero temporal apunta a NULL, se modifica "temp_storage_bytes" con el
//   tamaño de memoria temporal requerida */

//   cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
//       morton, morton_out, indices, indices_out, N);
//   /* Allocate temporary storage */
//   cudaMallocManaged(&d_temp_storage, temp_storage_bytes);
//   cudaMemPrefetchAsync(d_temp_storage, temp_storage_bytes, deviceId);
  
//   /* solo ordeno los índices */
//   cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
//       morton, morton_out, indices, indices_out, N);

//   mortonError = cudaGetLastError();
//   if(mortonError != cudaSuccess) printf("Error SORT: %s\n", cudaGetErrorString(mortonError));

//   asyncErr = cudaDeviceSynchronize();
//   if(asyncErr != cudaSuccess) printf("Error SORT sync: %s\n", cudaGetErrorString(asyncErr));

//   mytime = cast_t(tempo_t::now() - i_start).count();
//   std::cout << "  SORT time elapsed: " << mytime << " ms\n";
//   totaltime += mytime;

// #ifdef DEBUG
//   for(int i = 0; i<N; i++){
//     std::cout << morton_out[i] << " ";
//   }
//   std::cout << "\n";
//   for(int i = 0; i<N; i++){
//     std::cout << indices_out[i] << " ";
//   }
//   std::cout << "\n";
// #endif

//   i_start = tempo_t::now();

//   cudaMemPrefetchAsync(morton64, size_morton64, deviceId);

//   check_morton<<<numberOfBlocks, threadsPerBlock>>>( morton, indices, morton64, N );

//   mortonError = cudaGetLastError();
//   if(mortonError != cudaSuccess) printf("Error CHECK: %s\n", cudaGetErrorString(mortonError));

//   asyncErr = cudaDeviceSynchronize();
//   if(asyncErr != cudaSuccess) printf("Error CHECK sync: %s\n", cudaGetErrorString(asyncErr));

//   mytime = cast_t(tempo_t::now() - i_start).count();
//   std::cout << "  CHECK time elapsed: " << mytime << " ms\n";
//   totaltime += mytime;


//   // node_t default_node;
//   // default_node.parent_idx = 0xFFFFFFFF;
//   // default_node.left_idx   = 0xFFFFFFFF;
//   // default_node.right_idx  = 0xFFFFFFFF;
//   // default_node.object_idx = 0xFFFFFFFF;

//   // i_start = tempo_t::now();

//   // cudaMemset( node_list, 0xFFFFFFFF, size_nodes );

//   // mytime = cast_t(tempo_t::now() - i_start).count();
//   // std::cout << "  MemSet time elapsed: " << mytime << " ms\n";

//   i_start = tempo_t::now();

//   // cudaMemsetAsync( node_list, 0xFFFFFFFF, size_nodes*0.5 , (cudaStream_t)1);

//   cudaMemPrefetchAsync(aux_node_list, size_nodes*0.5, deviceId);
//   cudaMemPrefetchAsync(aux_aabb_list, size_aabbs*0.5, deviceId);

//   mortonError = cudaGetLastError();
//   if(mortonError != cudaSuccess) printf("Error LEAFS memory: %s\n", cudaGetErrorString(mortonError));

//   init_leafs<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, aux_aabb_list, point_cloud, indices, N );

//   // init_leafs<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, indices, N );

//   // init_leafs2<<<numberOfBlocks, threadsPerBlock>>>( node_list, indices, N, NN );

//   // init_indices<<<numberOfBlocks, threadsPerBlock>>>( aux_node_list, indices, N );

//   mortonError = cudaGetLastError();
//   if(mortonError != cudaSuccess) printf("Error LEAFS: %s\n", cudaGetErrorString(mortonError));

//   asyncErr = cudaDeviceSynchronize();
//   if(asyncErr != cudaSuccess) printf("Error LEAFS sync: %s\n", cudaGetErrorString(asyncErr));

//   mytime = cast_t(tempo_t::now() - i_start).count();
//   std::cout << "  LEAFS and AABBs time elapsed: " << mytime << " ms\n";
//   totaltime += mytime;

//   i_start = tempo_t::now();

//   cudaMemPrefetchAsync(node_list, size_nodes, deviceId);
//   cudaMemPrefetchAsync(morton64, size_morton64, deviceId);

//   // init_nodes<<<numberOfBlocks, threadsPerBlock>>>( node_list, morton64, N );

//   thrust::for_each(thrust::device,
//       thrust::make_counting_iterator<uint32_t>(0),
//       thrust::make_counting_iterator<uint32_t>(N - 1),
//       [node_list, morton64, N] __device__ (const uint32_t idx)
//       {
//         node_t aux;

//         aux.object_idx = 0xFFFFFFFF; //  internal nodes

//         const uint2 ij  = determine_range(morton64, N, idx);
//         const int gamma = find_split(morton64, ij.x, ij.y);

//         aux.left_idx  = gamma;
//         aux.right_idx = gamma + 1;
//         if(thrust::min(ij.x, ij.y) == gamma)
//         {
//             aux.left_idx += N - 1;
//         }
//         if(thrust::max(ij.x, ij.y) == gamma + 1)
//         {
//             aux.right_idx += N - 1;
//         }
//         node_list[aux.left_idx].parent_idx  = idx;
//         node_list[aux.right_idx].parent_idx = idx;

//         node_list[idx] = aux;
//         return;
//       });

//   mortonError = cudaGetLastError();
//   if(mortonError != cudaSuccess) printf("Error NODES: %s\n", cudaGetErrorString(mortonError));

//   // asyncErr = cudaDeviceSynchronize();
//   // if(asyncErr != cudaSuccess) printf("Error NODES sync: %s\n", cudaGetErrorString(asyncErr));

//   mytime = cast_t(tempo_t::now() - i_start).count();
//   std::cout << "  INTERNAL NODES time elapsed: " << mytime << " ms\n";
//   totaltime += mytime;


//   i_start = tempo_t::now();

//   cudaMemPrefetchAsync(flags, N*sizeof(int), deviceId);
//   cudaMemPrefetchAsync(aabb_list, size_nodes, deviceId);

//   create_aabb<<<numberOfBlocks, threadsPerBlock>>>( node_list, aabb_list, flags, N, NN );

//   // thrust::for_each(thrust::device,
//   //     thrust::make_counting_iterator<uint32_t>(N),
//   //     thrust::make_counting_iterator<uint32_t>(NN),
//   //     [node_list, aabb_list, flags] __device__ (const uint32_t idx)
//   //     {
//   //         uint32_t parent = node_list[idx].parent_idx;
//   //         while(parent != 0xFFFFFFFF) // means idx == 0
//   //         {
//   //             const int old = atomicCAS(flags + parent, 0, 1);
//   //             if(old == 0)
//   //             {
//   //                 // this is the first thread entered here.
//   //                 // wait the other thread from the other child node.
//   //                 return;
//   //             }
//   //             assert(old == 1);
//   //             // here, the flag has already been 1. it means that this
//   //             // thread is the 2nd thread. merge AABB of both childlen.

//   //             const uint32_t lidx = node_list[parent].left_idx;
//   //             const uint32_t ridx = node_list[parent].right_idx;
//   //             const aabb_t lbox = aabb_list[lidx];
//   //             const aabb_t rbox = aabb_list[ridx];
//   //             aabb_list[parent] = lbvh::merge(lbox, rbox);

//   //             // look the next parent...
//   //             parent = node_list[parent].parent_idx;
//   //         }
//   //         return;
//   //     });

//   mortonError = cudaGetLastError();
//   if(mortonError != cudaSuccess) printf("Error AABBs: %s\n", cudaGetErrorString(mortonError));

//   asyncErr = cudaDeviceSynchronize();
//   if(asyncErr != cudaSuccess) printf("Error AABBs sync: %s\n", cudaGetErrorString(asyncErr));

//   mytime = cast_t(tempo_t::now() - i_start).count();
//   std::cout << "  Create AABB time elapsed: " << mytime << " ms\n";
//   totaltime += mytime;
//   std::cout << "  CREATION takes: " << totaltime << " ms\n\n\n";


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


  return 0;
}
