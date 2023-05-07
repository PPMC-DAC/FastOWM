#include <thrust/sort.h>

#ifdef CHECK
#include <random>
#endif

namespace sycl_builder {

__global__ 
void init_node_list( node_t* nodes, aabb_t* bb, const uint32_t num_nodes,
        const node_t default_node, const aabb_t default_aabb)
{
  uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_nodes; i += stride)
  {

    nodes[i] = default_node;

    // n[i].parent_idx = 0xFFFFFFFFu;
    // n[i].left_idx = 0xFFFFFFFFu;
    // n[i].right_idx = 0xFFFFFFFFu;
    // n[i].object_idx = 0xFFFFFFFFu;
    // n[i].numPts = 0xFFFFFFFFu;
    // n[i].min_idx = 0xFFFFFFFFu;

    bb[i] = default_aabb;

  }
  return;
}

__global__ 
void init_morton_indices( uint32_t* morton, const uint32_t num_objects)
{
  uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_objects; i += stride)
  {

    morton[i] = 0u;

    // indices[i] = i;

  }
  return;
}

__global__ 
void init_flags( uint32_t* flags, const uint32_t num_internal_nodes)
{
  uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for(uint32_t i = index; i<num_internal_nodes; i += stride)
  {

    flags[i] = 0u;

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
void init_leafs_size(node_t* n, aabb_t* bb, const point_t* p,
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
void init_internal_nodes( node_t* n, const uint32_t num_leafs )
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



__global__
void create_aabb_size( node_t* n, aabb_t* bb, uint32_t* flags, const point_t* p, 
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




void SYCL_builder::reset()
{


  // cudaOccupancyMaxPotentialBlockSize( &numberOfBlocks,
  //                                     &threadsPerBlock,
  //                                     sycl_builder::init_struct );
  numberOfBlocks = 60;
  threadsPerBlock = 128;

	// cudaError_t lastError = cudaGetLastError();
	// if(lastError != cudaSuccess) printf("Error OCCUPANCY: %s\n", cudaGetErrorString(lastError));
	// cudaError_t syncError = cudaDeviceSynchronize();
	// if(syncError != cudaSuccess) printf("Sync error OCCUPANCY: %s\n", cudaGetErrorString(syncError));

  // sycl::host_accessor hm(morton);

  // for(int i=0; i<numObjects; i++){
  //   hm[i] = 0u;
  // }

  // sycl::host_accessor hf(flags);

  // for(int i=0; i<numInternalNodes; i++){
  //   hf[i] = 0u;
  // }


  device_queue.submit([&](sycl::handler& h) {
    auto accN = node_list.get_access<sycl::access::mode::write>(h);
    auto accA = aabb_list.get_access<sycl::access::mode::write>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dN = reinterpret_cast<node_t*>(ih.get_mem<sycl::backend::cuda>(accN));
      auto dA = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accA));

      sycl_builder::init_node_list<<<numberOfBlocks, threadsPerBlock>>>( 
                                    dN, dA, numNodes, default_node, default_aabb);
    });
  });

  device_queue.submit([&](sycl::handler& h) {
    auto accM = morton.get_access<sycl::access::mode::write>(h);
    // auto accI = indices.get_access<sycl::access::mode::write>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dM = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accM));
      // auto dI = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accI));

      sycl_builder::init_morton_indices<<<numberOfBlocks, threadsPerBlock>>>( dM, numObjects);
    });
  });

  device_queue.submit([&](sycl::handler& h) {
    auto accF = flags.get_access<sycl::access::mode::write>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dF = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accF));

      sycl_builder::init_flags<<<numberOfBlocks, threadsPerBlock>>>( dF, numInternalNodes);
    });
  });

	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error INIT: %s\n", cudaGetErrorString(lastError));
	cudaError_t syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error INIT: %s\n", cudaGetErrorString(syncError));


  // auto hn = node_list.get_access<sycl::access::mode::read>();

  // for(int i=numInternalNodes; i<numNodes; i++)
  // {
  //   if(hn[i].left_idx != 0xFFFFFFFFu || hn[i].right_idx != 0xFFFFFFFFu) printf("ERROR LEAF pos %d / %d : %u %u %u\n", i, numNodes, hn[i].left_idx, hn[i].right_idx, 0xFFFFFFFFu);
  // }

} // reset()




SYCL_builder::SYCL_builder( std::string inputTXT, const uint32_t chunk, cl::sycl::queue& q,
                    const uint32_t no, const uint32_t nn, const uint32_t nin) : 
            leafSize(chunk), device_queue(q), point_cloud(no), ord_point_cloud(no), morton(no),
            indices(no), node_list(nn), aabb_list(nn), flags(nin)
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

  /*accesor para obtener los puntos del fichero de entrada*/
  sycl::host_accessor hp(point_cloud);

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, NUM_PROCS);
  
#ifdef DEBUG
  std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif

#ifndef CHECK

  map_file_atof_tbb_th(inputTXT, hp, numObjects);

#else

  for(int i=0; i<numObjects; i++){
    hp[i].x = (int)uni(mt);
    hp[i].y = (int)uni(mt);
    hp[i].z = (int)uni(mt);
  }
  /*este tiene que ser siempre el mínimo en root si todo ha ido bien*/
  hp[9].z = -1.0;

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



  device_queue.submit([&](sycl::handler& h) {
    auto accP = point_cloud.get_access<sycl::access::mode::read>(h);
    auto accM = morton.get_access<sycl::access::mode::write>(h);
    auto accI = indices.get_access<sycl::access::mode::write>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dP = reinterpret_cast<point_t*>(ih.get_mem<sycl::backend::cuda>(accP));
      auto dM = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accM));
      auto dI = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accI));

      sycl_builder::get_morton_code<<<numberOfBlocks, threadsPerBlock>>>( 
                          dP, 
                          BBox, 
                          diffBox, 
                          dM,
                          dI,
                          numObjects );    
      });
  });


	
#ifdef DEBUG
	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error MORTON: %s\n", cudaGetErrorString(lastError));
	cudaError_t syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error MORTON: %s\n", cudaGetErrorString(syncError));

  double mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " MORTON time elapsed: " << mytime << " ms\n";
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


#ifdef DEBUG
  i_start = tempo_t::now();
#endif



  device_queue.submit([&](sycl::handler& h) {
    auto accM = morton.get_access<sycl::access::mode::read_write>(h);
    auto accI = indices.get_access<sycl::access::mode::read_write>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dM = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accM));
      auto dI = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accI));

      thrust::stable_sort_by_key(thrust::device, 
        dM, 
        dM + numObjects,
        dI,
        thrust::less<uint32_t>()
      );

    });
  });


	
#ifdef DEBUG
  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " SORT time elapsed: " << mytime << " ms\n";
#endif


#ifdef DEBUG
  i_start = tempo_t::now();
#endif



  device_queue.submit([&](sycl::handler& h) {
    auto accI = indices.get_access<sycl::access::mode::read>(h);
    auto accP = point_cloud.get_access<sycl::access::mode::read>(h);
    auto accO = ord_point_cloud.get_access<sycl::access::mode::write>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dI = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accI));
      auto dP = reinterpret_cast<point_t*>(ih.get_mem<sycl::backend::cuda>(accP));
      auto dO = reinterpret_cast<point_t*>(ih.get_mem<sycl::backend::cuda>(accO));

      sycl_builder::order_points<<<numberOfBlocks, threadsPerBlock>>>(
                        dP,
                        dO,
                        dI,
                        numObjects );
    });
  });



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



  device_queue.submit([&](sycl::handler& h) {
    auto accN = node_list.get_access<sycl::access::mode::write>(h);
    auto accA = aabb_list.get_access<sycl::access::mode::write>(h);
    auto accO = ord_point_cloud.get_access<sycl::access::mode::read>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dN = reinterpret_cast<node_t*>(ih.get_mem<sycl::backend::cuda>(accN));
      auto dA = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accA));
      auto dO = reinterpret_cast<point_t*>(ih.get_mem<sycl::backend::cuda>(accO));

      sycl_builder::init_leafs_size<<<numberOfBlocks, threadsPerBlock>>>( 
                        &dN[numInternalNodes], 
                        &dA[numInternalNodes], 
                        dO, 
                        numObjects,
                        numLeafs,
                        leafSize,
                        default_aabb);
    });
  });



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



  device_queue.submit([&](sycl::handler& h) {
    auto accN = node_list.get_access<sycl::access::mode::write>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dN = reinterpret_cast<node_t*>(ih.get_mem<sycl::backend::cuda>(accN));

      sycl_builder::init_internal_nodes<<<numberOfBlocks, threadsPerBlock>>>( 
                        dN, 
                        numLeafs );
    });
  });



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
									


  device_queue.submit([&](sycl::handler& h) {
    auto accN = node_list.get_access<sycl::access::mode::read_write>(h);
    auto accA = aabb_list.get_access<sycl::access::mode::read_write>(h);
    auto accF = flags.get_access<sycl::access::mode::read_write>(h);
    auto accO = ord_point_cloud.get_access<sycl::access::mode::read>(h);

    h.interop_task([=](sycl::interop_handler ih) {
      auto dN = reinterpret_cast<node_t*>(ih.get_mem<sycl::backend::cuda>(accN));
      auto dA = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accA));
      auto dF = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accF));
      auto dO = reinterpret_cast<point_t*>(ih.get_mem<sycl::backend::cuda>(accO));

      sycl_builder::create_aabb_size<<<numberOfBlocks, threadsPerBlock>>>( 
                        dN, 
                        dA, 
                        dF,
                        dO, 
                        numInternalNodes, 
                        numNodes );
    });
  });


	
#ifdef DEBUG
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error AABBs: %s\n", cudaGetErrorString(lastError));
	syncError = cudaDeviceSynchronize();
	if(syncError != cudaSuccess) printf("Sync error AABBs: %s\n", cudaGetErrorString(syncError));

  mytime = cast_t(tempo_t::now() - i_start).count();
  std::cout << " AABBs time elapsed: " << mytime << " ms\n";
#endif

  // auto hn = node_list.get_access<sycl::access::mode::read>();

  // for(int i=0; i<numInternalNodes; i++)
  // {
  //   if(hn[i].object_idx != 0xFFFFFFFFu) printf("ERROR INTER pos %d / %d : %u %u\n", i, numInternalNodes, hn[i].object_idx, 0xFFFFFFFFu);
  // }

  // for(int i=numInternalNodes; i<numNodes; i++)
  // {
  //   if(hn[i].left_idx != 0xFFFFFFFFu || hn[i].right_idx != 0xFFFFFFFFu) printf("ERROR LEAF pos %d / %d : %u %u %u\n", i, numNodes, hn[i].left_idx, hn[i].right_idx, 0xFFFFFFFFu);
  // }


#ifdef CHECK
  auto hn = node_list.get_access<sycl::access::mode::read>();
  auto ha = aabb_list.get_access<sycl::access::mode::read>();
  auto hp = ord_point_cloud.get_access<sycl::access::mode::read>();
  auto hf = flags.get_access<sycl::access::mode::read>();

  std::cout << std::endl;
  for(int i = 0; i<numNodes; i++){
    std::cout << hn[i].parent_idx << ", " << hn[i].left_idx << ", ";
    std::cout << hn[i].right_idx << ", " << hn[i].object_idx << ", ";
    std::cout << hn[i].numPts << ", " << hn[i].min_idx << "(";
    if(hn[i].min_idx != 0xFFFFFFFFu)
      std::cout << hp[hn[i].min_idx].z;
    std::cout << ")\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<numNodes; i++){
    std::cout << ha[i].upper.x << "," << ha[i].upper.y << " ";
    std::cout << ha[i].lower.x << "," << ha[i].lower.y << "\n";
  }
  std::cout << std::endl;
  for(int i = 0; i<numInternalNodes; i++){
    std::cout << hf[i] << " ";
  }
  std::cout << "\n\n";
#endif

  return;


} // build()


