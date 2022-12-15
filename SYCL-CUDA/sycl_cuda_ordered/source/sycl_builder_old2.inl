// #include <cub/cub.cuh>
// #include <cub/device/device_radix_sort.cuh>

// /*Esto evita un problema de compilación con thrust y el método
// va_printf. Consultar los comentarios en:
// https://github.com/NVIDIA/thrust/issues/1032*/
// using namespace cub;

// #include <thrust/pair.h>
// #include <thrust/tuple.h>
// #include <thrust/host_vector.h>
// #include <thrust/device_vector.h>
// #include <thrust/functional.h>
// #include <thrust/scan.h>
#include <thrust/sort.h>
// #include <thrust/fill.h>
// #include <thrust/for_each.h>
// #include <thrust/transform.h>
// #include <thrust/reduce.h>
// #include <thrust/iterator/constant_iterator.h>
// #include <thrust/iterator/counting_iterator.h>
// #include <thrust/execution_policy.h>

#ifdef CHECK
#include <random>
#endif

// using atomic_acc = sycl::ONEAPI::atomic_ref<uint32_t,   sycl::ONEAPI::memory_order::relaxed, 
//                                                         sycl::ONEAPI::memory_scope::system,
//                                                         sycl::access::address_space::global_space>;


namespace sycl_builder 
{

__global__ void init_struct( uint32_t* m, node_t* n, aabb_t* bb, uint32_t* f,
      const uint32_t num_objects, const uint32_t num_nodes, const uint32_t num_internal_nodes,
      const node_t default_node, const aabb_t default_aabb)
{
  uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_objects; i += stride)
  {
    m[i] = 0u;

    if(i >= num_nodes)
      return;

    n[i] = default_node;
    bb[i] = default_aabb;

    if(i >= num_internal_nodes)
      return;

    f[i] = 0u;
  }
  return;
}


__global__ 
void get_morton_code( const point_t* points, const whole_t BBox, const point_t diffBox, 
                        uint32_t* morton, uint32_t* indices, const uint32_t num_objects)
{
  uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_objects; i += stride){

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
void order_points(const point_t* p, point_t* op, const uint32_t* idx, const uint32_t num_objects)
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  const uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_objects; i += stride)
  {
    op[i] = p[idx[i]];
  }

  return;
}

/*con esta puedo evitar usar check_morton, porque a la misma vez inicializo los nodos hoja
voy calculando los aabb y el morton64 asociado*/
__global__
void init_leafs_size(bintree_node* n, aabb_t* bb, const point_t* p,
  const uint32_t num_objects, const uint32_t num_leafs, const uint32_t leafSize, const aabb_t default_aabb)
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  const uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_leafs; i += stride)
  {
  
    const uint32_t start = i*leafSize;
    const uint32_t end = i*leafSize+leafSize;
    // const uint32_t end = (start+leafSize <= num_objects)? start+leafSize : num_objects;
    aabb_t auxbb = default_aabb;

    real_t min = inf;
    uint32_t idmin = 0xFFFFFFFFu;

    /* modifico el índice del nodo hoja par que apunte al
    índice del primer punto en el vector ordenado de índices,
    es decir, este NO es el índice del primer punto sino el índice
    del primer índice en el vector de índices */
    n[i].object_idx = start;

    int j=start;
    for(; j<end && j<num_objects; j++) {
    // for(; j<end; j++) {
      // copy the point
      point_t auxp = p[j];
      // obtain the new aabb
      auxbb = merge(auxbb, {{auxp.x,auxp.y},{auxp.x,auxp.y}});
      // compare min
      if(auxp.z-min < TOL){
        min = auxp.z;
        idmin = j;
      }
    }

    __syncthreads();

    n[i].numPts = j-start;
    n[i].min_idx = idmin;
    bb[i] = auxbb;

    // __syncthreads();
  }

  return;
}


__global__
void init_nodes5( bintree_node* n, const uint32_t num_leafs )
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  const uint32_t stride = blockDim.x * gridDim.x;

  for(int i=index; i<num_leafs-1; i+=stride)
  {
    
    uint32_t idx = i;

    const uint64_t self_code = idx;
    const int L_delta = (i==0)? FW_S32_MIN : ::__clzll(self_code ^ (idx-1));
    const int R_delta = ::__clzll(self_code ^ (idx+1));
    const int d = (R_delta > L_delta) ? 1 : -1;

    // Compute upper bound for the length of the range

    const int delta_min = ::min(L_delta, R_delta);
    uint32_t l_max = 64;
    int i_tmp = idx + l_max * d;

    do{

      l_max<<=1;
      i_tmp = idx + l_max * d;

    } while( 0 <= i_tmp && i_tmp < num_leafs && delta_min < ::__clzll(self_code ^ i_tmp) );


    // Find the other end by binary search
    uint32_t l = 0;
    for(uint32_t t= l_max >> 1; t>0; t>>=1)
    {
      i_tmp = idx + (l + t) * d;
      if( 0 <= i_tmp && i_tmp < num_leafs && delta_min < ::__clzll(self_code ^ i_tmp))
        l += t;
    }

    uint32_t jdx = idx + l * d;
    if(d < 0)
    {
        swap(idx, jdx); // make it sure that idx < jdx
    }

    const uint64_t first_code = idx;
    const uint32_t prefix_node = ::__clzll(first_code ^ jdx);

    // binary search...
    uint32_t gamma  = idx;
    uint32_t stride = l;

    do
    {
        stride = (stride + 1) >> 1;
        const uint32_t middle = gamma + stride;
        if( middle < jdx && prefix_node < ::__clzll(first_code ^ middle))
          gamma = middle;

    } while(stride > 1);

    // k = gamma

    uint32_t lidx = (idx == gamma) ? (gamma + num_leafs - 1) : gamma;
    uint32_t ridx = (jdx == gamma + 1) ? (gamma + num_leafs) : (gamma + 1);

    n[i].left_idx = lidx;
    n[i].right_idx = ridx;
    n[lidx].parent_idx  = i;
    n[ridx].parent_idx = i;

    __syncthreads();
  }

  return;

}


/*
Igual que el caso del kernel anterior, en este tampoco puedo lanzar un thread por
cada nodo hoja y tengo que utilizar un stride para que cada thread haga más trabajo
*/
__global__
void create_aabb_size( bintree_node* n, aabb_t* bb, uint32_t* flags, const point_t* p, 
  const uint32_t num_internal_nodes, const uint32_t num_nodes)
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

    __syncthreads();
  }

  return;
}


void* mallocWrap(const size_t& size, cl::sycl::queue device_queue)
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


SYCL_builder::SYCL_builder( std::string inputTXT, const uint32_t chunk, cl::sycl::queue q) : 
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

  point_cloud = (point_t*)sycl_builder::mallocWrap(numObjects*sizeof(point_t), device_queue);
  ord_point_cloud = (point_t*)sycl_builder::mallocWrap(numObjects*sizeof(point_t), device_queue);

  morton = (uint32_t*)sycl_builder::mallocWrap(numObjects*sizeof(uint32_t), device_queue);
  indices = (uint32_t*)sycl_builder::mallocWrap(numObjects*sizeof(uint32_t), device_queue);

  node_list = (node_t*)sycl_builder::mallocWrap(numNodes * sizeof(node_t), device_queue);
  aabb_list = (aabb_t*)sycl_builder::mallocWrap(numNodes * sizeof(aabb_t), device_queue);

  flags = (uint32_t*)sycl_builder::mallocWrap(numInternalNodes*sizeof(uint32_t), device_queue);


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



void SYCL_builder::build()
{

#ifdef DEBUG
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif



	sycl_builder::get_morton_code<<<numberOfBlocks, threadsPerBlock>>>( 
                    point_cloud, 
                    BBox, 
                    diffBox, 
                    morton,
                    indices,
                    numObjects );

	
#ifdef DEBUG
	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error MORTON: %s\n", cudaGetErrorString(lastError));
	cudaError_t syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error MORTON: %s\n", cudaGetErrorString(syncError));

  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " MORTON time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif

// #ifdef CHECK
//   for(int i = 0; i<numObjects; i++){
//     std::cout << morton[i] << " ";
//   }
//   std::cout << "\n";
//   for(int i = 0; i<numObjects; i++){
//     std::cout << indices[i] << " ";
//   }
//   std::cout << "\n";
// #endif

  thrust::stable_sort_by_key(thrust::device, 
    morton, 
    morton + numObjects,
    indices,
    thrust::less<uint32_t>()
  );


	
#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " SORT time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif


  sycl_builder::order_points<<<numberOfBlocks, threadsPerBlock>>>(
                    point_cloud,
                    ord_point_cloud,
                    indices,
                    numObjects );


#ifdef DEBUG
  // e1.wait();
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error ORDER: %s\n", cudaGetErrorString(lastError));
	syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error ORDER: %s\n", cudaGetErrorString(syncError));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " ORDER time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif



	sycl_builder::init_leafs_size<<<numberOfBlocks, threadsPerBlock>>>( 
										&node_list[numInternalNodes], 
										&aabb_list[numInternalNodes], 
										ord_point_cloud, 
										numObjects,
										numLeafs,
										leafSize,
										default_aabb);


#ifdef DEBUG
  // e2.wait();
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error LEAFS: %s\n", cudaGetErrorString(lastError));
	syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error LEAFS: %s\n", cudaGetErrorString(syncError));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " LEAFS and AABBs time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif												  



	sycl_builder::init_nodes5<<<numberOfBlocks, threadsPerBlock>>>( 
										node_list, 
										numLeafs );


#ifdef DEBUG
  // e3.wait();
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error NODES: %s\n", cudaGetErrorString(lastError));
	syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error NODES: %s\n", cudaGetErrorString(syncError));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " INTERNAL NODES time elapsed: " << mytime << " ms\n";
#endif

// #ifdef CHECK
//   for(int i = 0; i<numNodes; i++){
//     std::cout << node_list[i].parent_idx << ", " << node_list[i].left_idx << ", ";
//     std::cout << node_list[i].right_idx << ", " << node_list[i].object_idx << ", ";
//     std::cout << node_list[i].numPts << ", " << node_list[i].min_idx << "\n";
//   }
//   std::cout << std::endl;
//   for(int i = 0; i<numNodes; i++){
//       std::cout << aabb_list[i].upper.x << "," << aabb_list[i].upper.y << " ";
//       std::cout << aabb_list[i].lower.x << "," << aabb_list[i].lower.y << "\n";
//     }
//   std::cout << std::endl;
//   for(int i = 0; i<numInternalNodes; i++){
//     std::cout << flags[i] << " ";
//   }
//   std::cout << "\n\n";
// #endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif


	sycl_builder::create_aabb_size<<<numberOfBlocks, threadsPerBlock>>>( 
										node_list, 
										aabb_list, 
										flags,
										ord_point_cloud, 
										numInternalNodes, 
										numNodes );

	
#ifdef DEBUG
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error AABBs: %s\n", cudaGetErrorString(lastError));
	syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error AABBs: %s\n", cudaGetErrorString(syncError));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " AABBs time elapsed: " << mytime << " ms\n";
#endif


#ifdef CHECK
  node_t* node_list_h = (node_t*)std::malloc(numNodes * sizeof(node_t));
  point_t* ord_point_cloud_h = (point_t*)std::malloc(numObjects * sizeof(point_t));
  aabb_t* aabb_list_h = (aabb_t*)std::malloc(numNodes * sizeof(aabb_t));
  uint32_t* flags_h = (uint32_t*)std::malloc(numInternalNodes*sizeof(uint32_t));

  device_queue.memcpy(node_list_h, node_list, numNodes*sizeof(node_t)).wait();
  device_queue.memcpy(ord_point_cloud_h, ord_point_cloud, numObjects*sizeof(point_t)).wait();
  device_queue.memcpy(aabb_list_h, aabb_list, numNodes*sizeof(aabb_t)).wait();
  device_queue.memcpy(flags_h, flags, numInternalNodes*sizeof(uint32_t)).wait();

  std::cout << std::endl;
  for(int i = 0; i<numNodes; i++){
    std::cout << node_list_h[i].parent_idx << ", " << node_list_h[i].left_idx << ", ";
    std::cout << node_list_h[i].right_idx << ", " << node_list_h[i].object_idx << ", ";
    std::cout << node_list_h[i].numPts << ", " << node_list_h[i].min_idx << "(";
    if(node_list_h[i].min_idx != 0xFFFFFFFFu)
      std::cout << ord_point_cloud_h[node_list_h[i].min_idx].z;
    std::cout << ")\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<numNodes; i++){
    std::cout << aabb_list_h[i].upper.x << "," << aabb_list_h[i].upper.y << " ";
    std::cout << aabb_list_h[i].lower.x << "," << aabb_list_h[i].lower.y << "\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<numInternalNodes; i++){
    std::cout << flags_h[i] << " ";
  }
  std::cout << "\n\n";

  free(node_list_h);
  free(ord_point_cloud_h);
  free(aabb_list_h);
  free(flags_h);
#endif

  return;


} // build()


void SYCL_builder::reset()
{

  cudaOccupancyMaxPotentialBlockSize( &numberOfBlocks,
                                      &threadsPerBlock,
                                      sycl_builder::init_struct );

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error OCCUPANCY: %s\n", cudaGetErrorString(lastError));
	cudaError_t syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error OCCUPANCY: %s\n", cudaGetErrorString(syncError));


  sycl_builder::init_struct<<<numberOfBlocks, threadsPerBlock>>>( 
                    morton,
                    node_list, 
                    aabb_list, 
                    flags,
                    numObjects,
                    numNodes,
                    numInternalNodes, 
                    default_node, 
                    default_aabb );


	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error INIT: %s\n", cudaGetErrorString(lastError));
	syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error INIT: %s\n", cudaGetErrorString(syncError));


} // reset()