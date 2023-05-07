#include <cub/cub.cuh>
#include <cub/device/device_radix_sort.cuh>


namespace bintree64_builder {


__global__ void init_struct( bintree_node* n, aabb_t* bb, const size_t NN, 
    const bintree_node default_node, const aabb_t default_aabb)
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

__global__ 
void get_morton_code(point_t* points, whole_t BBox, point_t diffBox, 
                        uint64_t* morton, size_t num_objects)
{
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for(int i = index; i<num_objects; i += stride){

    point_t aux = points[i];

    aux.x -= BBox.lower.x;
    aux.y -= BBox.lower.y;
    aux.z -= BBox.lower.z;
    aux.x /= diffBox.x;
    aux.y /= diffBox.y;
    aux.z /= diffBox.z;

    morton[i] = morton3D(aux);

  }

  return;
}

/*con esta puedo evitar usar check_morton, porque a la misma vez inicializo los nodos hoja
voy calculando los aabb y el morton64 asociado*/
__global__
void init_leafs_size(bintree_node* n, aabb_t* bb, uint64_t* m, uint64_t* m64,
  const int num_objects, const int num_leafs, const int leafSize, const aabb_t default_aabb,
  const whole_t BBox, const point_t diffBox)
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  const uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_leafs; i += stride)
  {
  
    uint32_t start = i*leafSize;
    uint32_t end = i*leafSize+leafSize;
    aabb_t auxbb = default_aabb;

    real_t min = inf;
    uint32_t idmin = 0xFFFFFFFFu;

    /* modifico el índice del nodo hoja par que apunte al
    índice del primer punto en el vector ordenado de códigos Morton64,
    es decir, este NO es el índice del primer punto sino el índice
    del primer código Morton, que hay que decodificar para obtener el punto */
    n[i].object_idx = start;

    int j=start;
    for(; j<end && j<num_objects; j++) {
      // copy the point
      point_t auxp = point3D(m[j]);

      auxp.x = auxp.x * diffBox.x + BBox.lower.x;
      auxp.y = auxp.y * diffBox.y + BBox.lower.y;
      auxp.z = auxp.z * diffBox.z + BBox.lower.z;

      // obtain the new aabb
      auxbb = merge(auxbb, {{auxp.x,auxp.y},{auxp.x,auxp.y}});
      // compare min
      if(auxp.z-min < TOL){
        min = auxp.z;
        idmin = j;
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

    // real_v c = getCenter(auxbb);

    // c.x -= BBox.lower.x;
    // c.y -= BBox.lower.y;
    // c.x /= diffBox.x;
    // c.y /= diffBox.y;

    // uint64_t aux64 = morton2D(c);
    // aux64 <<= 32;
    // aux64 |= i;
    // m64[i] = aux64;

    /*Esta línea es equivalente a todo lo que estaba haciendo desde getCenter()
    hasta aquí, porque no estaba dándole ninguna resolución a los códigos Morton*/
    m64[i] = uint64_t(i); 

    n[i].numPts = j-start;
    n[i].min_idx = idmin;
    bb[i] = auxbb;

    __syncthreads();
  }

  return;
}

/*
Tengo que utilizar este kernel porque parece ser que hay un problema con el número
máximo de threads que puedo tener en vuelo, y no puedo lanzar un thread por cada nodo interno.
Los errores son de este tipo:
"CUDA error 700 [/usr/local/cuda/bin/../targets/x86_64-linux/include/cub/device/dispatch/
../../agent/../block/../iterator/../util_device.cuh, 630]: an illegal memory access was encountered"
*/
__global__
void init_nodes( bintree_node* n, const uint64_t* node_code, const uint32_t num_leafs )
{
  const uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  const uint32_t stride = blockDim.x * gridDim.x;

  for(int i=index; i<num_leafs-1; i+=stride)
  {
    
    uint32_t idx = i;

    const uint64_t self_code = node_code[idx];
    const int L_delta = (i==0)? FW_S32_MIN : ::__clzll(self_code ^ node_code[idx-1]);
    const int R_delta = ::__clzll(self_code ^ node_code[idx+1]);
    const int d = (R_delta > L_delta) ? 1 : -1;

    // Compute upper bound for the length of the range

    const int delta_min = ::min(L_delta, R_delta);
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
void create_aabb_size( bintree_node* n, aabb_t* bb, uint32_t* flags, const uint64_t* m, 
  const uint32_t num_internal_nodes, const uint32_t num_nodes, const whole_t BBox, const point_t diffBox)
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

      real_t lz = point3D(m[lmin]).z;
      real_t rz = point3D(m[rmin]).z;

      lz = lz * diffBox.z + BBox.lower.z;
      rz = rz * diffBox.z + BBox.lower.z;
      
      // compare mins
      if(lz - rz < TOL)
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


} // namespace bintree64_builder


Bintree64_builder::Bintree64_builder( std::string inputTXT, const uint32_t chunk) : leafSize(chunk)
{
  readHeader(inputTXT, BBox, numObjects);

  numLeafs = (uint32_t)((numObjects-1)/leafSize) + 1u;
  numInternalNodes = numLeafs - 1u;
  numNodes = 2*numLeafs - 1u;
  
  diffBox.x = BBox.upper.x - BBox.lower.x;
  diffBox.y = BBox.upper.y - BBox.lower.y;
  diffBox.z = BBox.upper.z - BBox.lower.z;

  cudaGetDevice(&deviceId);
  cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);

  cudaMallocManaged((void**)&point_cloud, numObjects*sizeof(point_t));
  cudaMallocManaged((void**)&morton, numObjects*sizeof(uint64_t));
  cudaMallocManaged((void**)&morton_out, numObjects*sizeof(uint64_t));
  cudaMallocManaged((void**)&morton64, numLeafs*sizeof(uint64_t));
  cudaMallocManaged((void**)&node_list, numNodes*sizeof(bintree_node));
  cudaMallocManaged((void**)&aabb_list, numNodes*sizeof(aabb_t));
  cudaMallocManaged((void**)&flags, numInternalNodes*sizeof(uint32_t));

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
                                      bintree64_builder::get_morton_code );

  bintree64_builder::init_struct<<<numberOfBlocks, threadsPerBlock>>>( 
                    node_list, 
                    aabb_list, 
                    numNodes, 
                    default_node, 
                    default_aabb );

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error INIT NODES and ABBBs: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

}


void Bintree64_builder::prefetch()
{
  return;
}

template <typename pointer, typename size, typename ...Rest>
void Bintree64_builder::prefetch( pointer p, size s, int device, Rest... args )
{
  cudaMemPrefetchAsync(p, s, device);
  prefetch( args... );
  return;
}


void Bintree64_builder::prefetchAll()
{
	cudaMemPrefetchAsync((void*)node_list,   numNodes*sizeof(bintree_node),     deviceId);
	cudaMemPrefetchAsync((void*)aabb_list,   numNodes*sizeof(aabb_t),           deviceId);
	cudaMemPrefetchAsync((void*)point_cloud, numObjects*sizeof(point_t),        deviceId);
	cudaMemPrefetchAsync((void*)morton,      numObjects*sizeof(uint64_t),       deviceId);
	cudaMemPrefetchAsync((void*)morton_out,  numObjects*sizeof(uint64_t),       deviceId);
	cudaMemPrefetchAsync((void*)morton64,    numLeafs*sizeof(uint64_t),         deviceId);
	cudaMemPrefetchAsync(flags,       numInternalNodes*sizeof(uint32_t), deviceId);

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error PREFETCH: %s\n", cudaGetErrorString(lastError));

  return;
}


void Bintree64_builder::build()
{

#ifdef DEBUG
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif

	bintree64_builder::get_morton_code<<<numberOfBlocks, threadsPerBlock>>>( 
                    point_cloud, 
                    BBox, 
                    diffBox, 
                    morton, 
                    numObjects );
	
	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error MORTON: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

#ifdef DEBUG
  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " MORTON time elapsed: " << mytime << " ms\n";

  i_start = tempo_t::now();
#endif


  size_t temp_storage_bytes = 0;
	/* Determine temporary device storage requirements */
	cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes,
		morton, morton_out, numObjects);

	/* Allocate temporary device storage */
	cudaMalloc(&d_temp_storage, temp_storage_bytes);
	
	/* Sort indices */
	cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes,
		morton, morton_out, numObjects);
	
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error SORT: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();


#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " SORT time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif

	bintree64_builder::init_leafs_size<<<numberOfBlocks, threadsPerBlock>>>( 
										&node_list[numInternalNodes], 
										&aabb_list[numInternalNodes], 
										morton_out,
										morton64, 
										numObjects,
										numLeafs,
										leafSize,
										default_aabb,
										BBox,
										diffBox );

	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error LEAFS: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " LEAFS and AABBs time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif
																		  
	bintree64_builder::init_nodes<<<numberOfBlocks, threadsPerBlock>>>( 
										node_list, 
										morton64, 
										numLeafs );

	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error NODES: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " INTERNAL NODES time elapsed: " << mytime << " ms\n";
#endif

#ifdef DEBUG
  i_start = tempo_t::now();
#endif
									
	bintree64_builder::create_aabb_size<<<numberOfBlocks, threadsPerBlock>>>( 
										node_list, 
										aabb_list, 
										flags,
										morton_out, 
										numInternalNodes, 
										numNodes,
                    BBox,
                    diffBox );
	
	
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error AABBs: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " AABBs time elapsed: " << mytime << " ms\n";
#endif

  return;
}

void Bintree64_builder::reset()
{
  cudaFree(d_temp_storage);
  d_temp_storage = NULL;

  cudaMemset(morton, 0u, numObjects*sizeof(uint64_t));
  cudaMemset(morton_out, 0u, numObjects*sizeof(uint64_t));
  // cudaMemset(indices, 0u, numObjects*sizeof(uint32_t));
  // cudaMemset(indices_out, 0u, numObjects*sizeof(uint32_t));
  cudaMemset(morton64, 0u, numLeafs*sizeof(uint64_t));
  cudaMemset(flags, 0u, numInternalNodes*sizeof(uint32_t));

  prefetchAll();

  cudaOccupancyMaxPotentialBlockSize( &numberOfBlocks,
                                      &threadsPerBlock,
                                      bintree64_builder::get_morton_code );

  bintree64_builder::init_struct<<<numberOfBlocks, threadsPerBlock>>>( 
                    node_list, 
                    aabb_list, 
                    numNodes, 
                    default_node, 
                    default_aabb );

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error RESET: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();
  
}

template<typename kernel_t>
int Bintree64_builder::getOccupancy( kernel_t kernel )
{
  int block_size; // max_blocks per SM

  cudaOccupancyMaxActiveBlocksPerMultiprocessor( &block_size, 
                                                 kernel,
                                                 threadsPerBlock, 
                                                 0 );

  return block_size * numberOfSMs;
}
