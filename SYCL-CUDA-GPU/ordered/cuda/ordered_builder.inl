#include <cub/cub.cuh>
#include <cub/device/device_radix_sort.cuh>


namespace ordered_builder {

__global__ void init_struct( bintree_node* n, aabb_t* bb, const uint32_t num_objects, 
    const bintree_node default_node, const aabb_t default_aabb)
{
  uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_objects; i += stride)
  {
    n[i] = default_node;
    bb[i] = default_aabb;
  }
  return;
}

__global__ void init_cloud( const point_t* cloud, const uint32_t* ids, point_t* point_cloud, const uint32_t num_objects)
{
  uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_objects; i += stride)
  {
    point_cloud[i] = cloud[ids[i]];
  }
  return;
}

__global__ 
void get_morton_code( const point_t* points, const whole_t BBox, const point_t diffBox, 
                        uint32_t* morton, uint32_t* indices, const uint32_t num_objects)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for(int i = index; i<num_objects; i += stride){

    indices[i] = i;

    point_t aux = points[i];

    __syncthreads();

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
    point_t auxp = p[idx[i]];
    __syncthreads();
    op[i] = auxp;
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


} // namespace ordered_builder


Ordered_builder::Ordered_builder( std::string inputTXT, const uint32_t chunk) : leafSize(chunk)
{
  readHeader(inputTXT, BBox, numObjects);

  numLeafs = (uint32_t)((numObjects-1)/leafSize) + 1u;
  numInternalNodes = numLeafs - 1u;
  numNodes = 2*numLeafs - 1u;
  
  diffBox.x = BBox.upper.x - BBox.lower.x;
  diffBox.y = BBox.upper.y - BBox.lower.y;

  cudaGetDevice(&deviceId);
  cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);

  cudaMallocManaged(&point_cloud, numObjects*sizeof(point_t));
  cudaMallocManaged(&ord_point_cloud, numObjects*sizeof(point_t));
  cudaMallocManaged(&morton, numObjects*sizeof(uint32_t));
  cudaMallocManaged(&morton_out, numObjects*sizeof(uint32_t));
  cudaMallocManaged(&indices, numObjects*sizeof(uint32_t));
  cudaMallocManaged(&indices_out, numObjects*sizeof(uint32_t));
  // cudaMallocManaged(&morton64, numLeafs*sizeof(uint64_t));
  cudaMallocManaged(&node_list, numNodes*sizeof(bintree_node));
  cudaMallocManaged(&aabb_list, numNodes*sizeof(aabb_t));
  cudaMallocManaged(&flags, numInternalNodes*sizeof(uint32_t));

  d_temp_storage = NULL;

  prefetch(point_cloud, numObjects*sizeof(point_t), cudaCpuDeviceId);

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, NUM_PROCS);
  
#ifdef DEBUG
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif
  map_file_atof_tbb_th(inputTXT, point_cloud, numObjects);
#ifdef DEBUG
  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << "Load time elapsed: " << mytime << " ms\n";
#endif

	prefetchAll();

  cudaOccupancyMaxPotentialBlockSize( &numberOfBlocks,
                                      &threadsPerBlock,
                                      ordered_builder::init_struct );

  // numberOfBlocks = 60;
  // threadsPerBlock = 128;

  ordered_builder::init_struct<<<numberOfBlocks, threadsPerBlock>>>( 
                    node_list, 
                    aabb_list, 
                    numNodes, 
                    default_node, 
                    default_aabb );

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error INIT NODES and ABBBs: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

}

Ordered_builder::Ordered_builder(const point_t* cloud, whole_t& BB, const uint32_t* minIDs, const uint32_t NObj, const uint32_t chunk) : 
                                  BBox(BB), numObjects(NObj), leafSize(chunk)
{

  numLeafs = (uint32_t)((numObjects-1)/leafSize) + 1u;
  numInternalNodes = numLeafs - 1u;
  numNodes = 2*numLeafs - 1u;
  
  diffBox.x = BBox.upper.x - BBox.lower.x;
  diffBox.y = BBox.upper.y - BBox.lower.y;

  cudaGetDevice(&deviceId);
  cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);

  cudaMallocManaged(&point_cloud, numObjects*sizeof(point_t));
  cudaMallocManaged(&ord_point_cloud, numObjects*sizeof(point_t));
  cudaMallocManaged(&morton, numObjects*sizeof(uint32_t));
  cudaMallocManaged(&morton_out, numObjects*sizeof(uint32_t));
  cudaMallocManaged(&indices, numObjects*sizeof(uint32_t));
  cudaMallocManaged(&indices_out, numObjects*sizeof(uint32_t));
  // cudaMallocManaged(&morton64, numLeafs*sizeof(uint64_t));
  cudaMallocManaged(&node_list, numNodes*sizeof(bintree_node));
  cudaMallocManaged(&aabb_list, numNodes*sizeof(aabb_t));
  cudaMallocManaged(&flags, numInternalNodes*sizeof(uint32_t));

  d_temp_storage = NULL;

	prefetchAll();

  // cudaMemPrefetchAsync(cloud, Npoints*sizeof(point_t), deviceId);
  // cudaMemPrefetchAsync(minIDs, numObjects*sizeof(uint32_t), deviceId);

  cudaOccupancyMaxPotentialBlockSize( &numberOfBlocks,
                                      &threadsPerBlock,
                                      ordered_builder::init_struct );

  // numberOfBlocks = 60;
  // threadsPerBlock = 128;

  ordered_builder::init_struct<<<numberOfBlocks, threadsPerBlock>>>( 
                    node_list, 
                    aabb_list, 
                    numNodes, 
                    default_node, 
                    default_aabb );

  ordered_builder::init_cloud<<<numberOfBlocks, threadsPerBlock>>>(
                    cloud,
                    minIDs,
                    point_cloud,
                    numObjects );

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error INIT_CLOUD: %s\n", cudaGetErrorString(lastError));
	lastError = cudaDeviceSynchronize();
	if(lastError != cudaSuccess) printf("Error INIT_CLOUD SYNC: %s\n", cudaGetErrorString(lastError));

}


void Ordered_builder::prefetch()
{
  return;
}

template <typename pointer, typename size, typename ...Rest>
void Ordered_builder::prefetch( pointer p, size s, int device, Rest... args )
{
  cudaMemPrefetchAsync(p, s, device);
  prefetch( args... );
  return;
}


void Ordered_builder::prefetchAll()
{
	cudaMemPrefetchAsync(point_cloud, numObjects*sizeof(point_t),        deviceId);
	cudaMemPrefetchAsync(ord_point_cloud, numObjects*sizeof(point_t),        deviceId);
	cudaMemPrefetchAsync(node_list,   numNodes*sizeof(bintree_node),     deviceId);
	cudaMemPrefetchAsync(aabb_list,   numNodes*sizeof(aabb_t),           deviceId);
	cudaMemPrefetchAsync(morton,      numObjects*sizeof(uint32_t),       deviceId);
	cudaMemPrefetchAsync(indices,     numObjects*sizeof(uint32_t),       deviceId);
	cudaMemPrefetchAsync(morton_out,  numObjects*sizeof(uint32_t),       deviceId);
	cudaMemPrefetchAsync(indices_out, numObjects*sizeof(uint32_t),       deviceId);
	// cudaMemPrefetchAsync(morton64,    numLeafs*sizeof(uint64_t),         deviceId);
	cudaMemPrefetchAsync(flags,       numInternalNodes*sizeof(uint32_t), deviceId);

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error PREFETCH: %s\n", cudaGetErrorString(lastError));

  return;
}


void Ordered_builder::build()
{

#ifdef DEBUG
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif

	ordered_builder::get_morton_code<<<numberOfBlocks, threadsPerBlock>>>( 
                    point_cloud, 
                    BBox, 
                    diffBox, 
                    morton, 
                    indices, 
                    numObjects );
	

#ifdef DEBUG
	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error MORTON: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " MORTON time elapsed: " << mytime << " ms\n";

  i_start = tempo_t::now();
#endif


  size_t temp_storage_bytes = 0;
	/* Determine temporary device storage requirements */
	cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
		morton, morton_out, indices, indices_out, numObjects);

	/* Allocate temporary device storage */
	cudaMalloc(&d_temp_storage, temp_storage_bytes);
	
	/* Sort indices */
	cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
		morton, morton_out, indices, indices_out, numObjects);

#ifdef DEBUG
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error SORT: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " SORT time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif

  ordered_builder::order_points<<<numberOfBlocks, threadsPerBlock>>>(
                    point_cloud,
                    ord_point_cloud,
                    indices_out,
                    numObjects );

#ifdef DEBUG
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error SORT: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();
  
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " ORDER time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif

	ordered_builder::init_leafs_size<<<numberOfBlocks, threadsPerBlock>>>( 
										&node_list[numInternalNodes], 
										&aabb_list[numInternalNodes], 
										ord_point_cloud, 
										// morton64, 
										numObjects,
										numLeafs,
										leafSize,
										default_aabb);


  /*En el momento que elimino la depenencia con morton64 estos dos kernels pueden ser concurrentes*/

#ifdef DEBUG
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error LEAFS: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " LEAFS and AABBs time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif
																		  
	ordered_builder::init_nodes5<<<numberOfBlocks, threadsPerBlock>>>( 
										node_list, 
										// morton64, 
										numLeafs );

#ifdef DEBUG
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error NODES: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " INTERNAL NODES time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif
									
	ordered_builder::create_aabb_size<<<numberOfBlocks, threadsPerBlock>>>( 
										node_list, 
										aabb_list, 
										flags,
										ord_point_cloud, 
										numInternalNodes, 
										numNodes );
	
#ifdef DEBUG
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error AABBs: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " AABBs time elapsed: " << mytime << " ms\n";
#endif

  return;
}

void Ordered_builder::reset()
{

  cudaFree(d_temp_storage);
  d_temp_storage = NULL;

  cudaMemset(morton, 0u, numObjects*sizeof(uint32_t));
  cudaMemset(morton_out, 0u, numObjects*sizeof(uint32_t));
  cudaMemset(indices, 0u, numObjects*sizeof(uint32_t));
  cudaMemset(indices_out, 0u, numObjects*sizeof(uint32_t));
  // cudaMemset(morton64, 0u, numLeafs*sizeof(uint64_t));
  cudaMemset(flags, 0u, numInternalNodes*sizeof(uint32_t));

  prefetchAll();

  // cudaOccupancyMaxPotentialBlockSize( &numberOfBlocks,
  //                                     &threadsPerBlock,
  //                                     ordered_builder::init_struct );

  ordered_builder::init_struct<<<numberOfBlocks, threadsPerBlock>>>( 
                    node_list, 
                    aabb_list, 
                    numNodes, 
                    default_node, 
                    default_aabb );

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error RESET: %s\n", cudaGetErrorString(lastError));
	cudaError_t syncError = cudaDeviceSynchronize();
	if(lastError != cudaSuccess) printf("Sync error RESET: %s\n", cudaGetErrorString(lastError));
  
}

template<typename kernel_t>
int Ordered_builder::getOccupancy( kernel_t kernel )
{
  int block_size; // max_blocks per SM

  cudaOccupancyMaxActiveBlocksPerMultiprocessor( &block_size, 
                                                 kernel,
                                                 threadsPerBlock, 
                                                 0 );

  return block_size * numberOfSMs;
}
