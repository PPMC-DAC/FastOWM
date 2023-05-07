#include <iostream>

#include <bintree/cuda/bintree_builder.h>

void bintree_test(std::string inputTXT, const uint32_t chunkDim, float& dtime)
{
	aabb_t BBox;
	uint32_t N;

	readHeader(inputTXT, BBox, N);
  
	real_v diffBox;
	diffBox.x = BBox.upper.x - BBox.lower.x;
	diffBox.y = BBox.upper.y - BBox.lower.y;  
  
	point_t*   point_cloud = NULL;
	uint32_t*  morton = NULL;
	uint32_t*  morton_out = NULL;
	uint32_t*  indices = NULL;
	uint32_t*  indices_out = NULL;
	void*      d_temp_storage = NULL;
	size_t temp_storage_bytes = 0;
	  
	uint64_t* morton64S = NULL;
	bintree_node* node_listS = NULL;
	aabb_t* aabb_listS = NULL;
	uint32_t* flagsS = NULL;
    
	int deviceId;
	int numberOfSMs;

	cudaGetDevice(&deviceId);
	cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId);
	// printf("Device ID: %d\tNumber of SMs: %d\n", deviceId, numberOfSMs);
	
	const uint32_t num_objects = N;
	/* reduced number of objects */
	const uint32_t r_num_objects = (uint32_t)((num_objects-1)/chunkDim) + 1u;
	const uint32_t r_num_internal_nodes = r_num_objects - 1u;
	const uint32_t r_num_nodes = 2*r_num_objects - 1u;

	cudaError_t lastError;
	// cudaError_t asyncErr;
  
	const uint32_t threadsPerBlock = 256;
	const uint32_t numberOfBlocks = 32*numberOfSMs;

	const uint32_t nBlocks_node_func = (uint32_t)((r_num_internal_nodes-1)/threadsPerBlock) + 1;
	const uint32_t nBlocks_aabb_func = (uint32_t)((r_num_objects-1)/threadsPerBlock) + 1;
  
	cudaMallocManaged(&point_cloud, num_objects*sizeof(point_t));
	cudaMallocManaged(&morton, num_objects*sizeof(uint32_t));
	cudaMallocManaged(&morton_out, num_objects*sizeof(uint32_t));
	cudaMallocManaged(&indices, num_objects*sizeof(uint32_t));
	cudaMallocManaged(&indices_out, num_objects*sizeof(uint32_t));
  
	cudaMallocManaged(&morton64S, r_num_objects*sizeof(uint64_t));
	cudaMallocManaged(&node_listS, r_num_nodes*sizeof(bintree_node));
	cudaMallocManaged(&aabb_listS, r_num_nodes*sizeof(aabb_t));
	cudaMallocManaged(&flagsS, r_num_internal_nodes*sizeof(uint32_t));

	cudaMemPrefetchAsync(node_listS, r_num_nodes*sizeof(bintree_node), deviceId);
	cudaMemPrefetchAsync(aabb_listS, r_num_nodes*sizeof(aabb_t), deviceId);
	cudaMemPrefetchAsync(morton, num_objects*sizeof(uint32_t), deviceId);
	cudaMemPrefetchAsync(indices, num_objects*sizeof(uint32_t), deviceId);
	cudaMemPrefetchAsync(morton_out, num_objects*sizeof(uint32_t), deviceId);
	cudaMemPrefetchAsync(indices_out, num_objects*sizeof(uint32_t), deviceId);
	cudaMemPrefetchAsync(morton64S, r_num_objects*sizeof(uint64_t), deviceId);
	cudaMemPrefetchAsync(flagsS, r_num_internal_nodes*sizeof(uint32_t), deviceId);
  
	cudaMemPrefetchAsync(point_cloud, num_objects*sizeof(point_t), cudaCpuDeviceId);

	if(read_pointsC(inputTXT, point_cloud, num_objects) < 0){
		printf("Unable to read file!\n");
		exit(-1);
	}

	cudaMemPrefetchAsync(point_cloud, num_objects*sizeof(point_t), deviceId);
  
	bintree_builder::init_struct<<<numberOfBlocks, threadsPerBlock>>>( node_listS, aabb_listS, r_num_nodes, default_node, default_aabb);

	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error INIT NODES y ABBBs: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

    cudaEvent_t start, stop;
    cudaEventCreate( &start );
    cudaEventCreate( &stop );

	cudaEventRecord( start, 0 );
	  
	bintree_builder::get_morton_code<<<numberOfBlocks, threadsPerBlock>>>( point_cloud, BBox, diffBox, morton, indices, num_objects);
	
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error MORTON: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

	/* Determine temporary device storage requirements */
	cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
		morton, morton_out, indices, indices_out, N);

	/* Allocate temporary storage */
	cudaMallocManaged(&d_temp_storage, temp_storage_bytes);
	cudaMemPrefetchAsync(d_temp_storage, temp_storage_bytes, deviceId);
	
	/* Sort indices */
	cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
		morton, morton_out, indices, indices_out, N);
	
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error SORT: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();

	bintree_builder::init_leafs_size2<<<numberOfBlocks, threadsPerBlock>>>( 
										&node_listS[r_num_internal_nodes], 
										&aabb_listS[r_num_internal_nodes], 
										point_cloud, 
										indices_out,
										morton64S, 
										num_objects,
										r_num_objects,
										chunkDim,
										default_aabb,
										BBox,
										diffBox );

	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error LEAFS: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();
																		  
	bintree_builder::init_nodes2<<<nBlocks_node_func, threadsPerBlock>>>( 
										node_listS, 
										morton64S, 
										r_num_objects );

	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error NODES: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();
									
	bintree_builder::create_aabb_size<<<nBlocks_aabb_func, threadsPerBlock>>>( 
										node_listS, 
										aabb_listS, 
										flagsS,
										point_cloud, 
										r_num_internal_nodes, 
										r_num_nodes );
	
	
	lastError = cudaGetLastError();
	if(lastError != cudaSuccess) printf("Error AABBs: %s\n", cudaGetErrorString(lastError));
	cudaDeviceSynchronize();
							
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );
	cudaEventElapsedTime( &dtime, start, stop );
	cudaEventDestroy( start );
	cudaEventDestroy( stop );

	cudaFree(point_cloud);
	cudaFree(morton);
	cudaFree(morton_out);
	cudaFree(indices);
	cudaFree(indices_out);
	cudaFree(d_temp_storage);
  
	cudaFree(morton64S);
	cudaFree(node_listS);
	cudaFree(aabb_listS);
	cudaFree(flagsS);
	  
	return;
}

int main(int argc, char* argv[])
{

	const std::string inputTXT = (argc > 1)? argv[1] : "../data/INAER_2011_Alcoy.xyz";
	const int chunkDim = (argc > 2)? atoi(argv[2]) : 256;

	float time = 0.0;
	// int nreps = 10;
	// for(int rep=0; rep<nreps; rep++)
	// {
		float dtime=0.0;
		bintree_test(inputTXT, chunkDim, dtime);
		time += dtime;
	// }

	std::cout << " CREATION takes: " << time << " ms\n";

	return 0;
}