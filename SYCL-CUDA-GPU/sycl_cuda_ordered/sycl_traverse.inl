#include <iostream>

#include <sycl_cuda_ordered/source/sycl_builder.h>

#include <basic/traverse_sycl_cuda.h>

// class CUDASelector : public sycl::device_selector {
// public:
//   int operator()(const sycl::device &Device) const override {
//     using namespace sycl::info;

//     const std::string DriverVersion = Device.get_info<device::driver_version>();

//     if (Device.is_gpu() && (DriverVersion.find("CUDA") != std::string::npos)) {
//       std::cout << "CUDA device found: ";
//       return 1;
//     };
//     return -1;
//   }
// };

void sycl_traverse(std::string inputTXT, uint32_t chunkDim)
{

#ifdef CHECK
    inputTXT = "CHECK";
    chunkDim = 2;
#endif

// #ifdef FEMU
//     cl::sycl::INTEL::fpga_emulator_selector device_selector{};
// #elif GPU
//     cl::sycl::gpu_selector device_selector{};
// #elif CPU
//     cl::sycl::cpu_selector device_selector{};
// #endif

    cl::sycl::queue device_queue(CUDASelector().select_device());

    std::cout << device_queue.get_device().get_info<sycl::info::device::name>()  
    << " @ " << device_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << "Mhz (" << 
    device_queue.get_device().get_info<sycl::info::device::max_compute_units>() << " cores)" << std::endl;

    uint32_t numObjects = 0;
    whole_t BBox;

#ifndef CHECK
    readHeader(inputTXT, BBox, numObjects);
#else
    numObjects = 32;
#endif

    uint32_t numLeafs = (uint32_t)((numObjects-1)/chunkDim) + 1u;
    uint32_t numInternalNodes = numLeafs - 1u;
    uint32_t numNodes = 2*numLeafs - 1u;

    SYCL_builder builder(inputTXT, chunkDim, device_queue, numObjects, numNodes, numInternalNodes);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.numObjects << " points; chunkDim: " << builder.leafSize << std::endl;
    std::cout << builder.numLeafs << " leaf nodes\n";
    std::cout << builder.numInternalNodes << " internal nodes\n";
    std::cout << builder.numNodes << " nodes\n\n";

    std::cout << "BBox: " << builder.BBox.upper.x << "," << builder.BBox.upper.y << " ";
    std::cout << builder.BBox.lower.x << "," << builder.BBox.lower.y << "\n\n";
  
    std::cout << "diffBox: " << builder.diffBox.x << "," << builder.diffBox.y << "\n\n";
#endif

    std::chrono::time_point<tempo_t> start = tempo_t::now();
    
    builder.build();

    cudaError_t lastError = cudaGetLastError();
    if(lastError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(lastError));
    cudaError_t syncError = cudaDeviceSynchronize();
    if(syncError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(syncError));

    std::double_t dtime = cast_t(tempo_t::now() - start).count();    
    std::cout << "  CREATION takes: " << dtime << " ms\n";


    uint32_t Wsize = 10;
    // uint32_t Bsize = 20;
    std::double_t Overlap = 0.8;

    uint32_t Ncells;
    uint32_t nRows, nCols;
    std::double_t Width, High, Density;

    Width = builder.diffBox.x;
    High = builder.diffBox.y;
    // Densidad en puntos/m^2
    Density = builder.numObjects/(Width*High);

    // El numero minimo sera la mitad del numero de puntos medio por celda
    uint32_t minNumPoints = (uint32_t)(0.5*Density*Wsize*Wsize);

    std::double_t Displace = round2d(Wsize*(1-Overlap));

    // Stage 1 parameters
    if(Overlap > 0.0) {
        nCols=(int)(round((Width+2*Wsize*Overlap)/Displace))-1;
        nRows=(int)(round((High+2*Wsize*Overlap)/Displace))-1;
    } else {
        nCols=(int)floor(Width/Wsize)+1;
        nRows=(int)floor(High/Wsize)+1;
    }
    Ncells = nCols*nRows;
#ifdef DEBUG
    printf("nCols: %d, nRows: %d\n", nCols, nRows);
    printf("minNumPoints: %d\n", minNumPoints);
#endif

    uint32_t* count = (uint32_t*)sycl_builder::mallocWrap(Ncells*sizeof(uint32_t), device_queue);

#if DEVICE
    uint32_t* count_h = (uint32_t*)std::malloc(Ncells*sizeof(uint32_t));
#endif

#ifndef DEBUG
    int n_tests = 40;
#else
    int n_tests = 1;
#endif

    std::double_t total = 0.0;

    for(int i=0; i<n_tests; i++){

#if DEVICE
        device_queue.memset(count, 0u, Ncells*sizeof(uint32_t)).wait();
#else
        std::memset(count, 0u, Ncells*sizeof(uint32_t));
#endif

        start = tempo_t::now();

        // stage1query(builder.node_list, builder.aabb_list, builder.ord_point_cloud, count,
        //         Wsize, Overlap, nCols, nRows, minNumPoints, builder.BBox, builder.numInternalNodes);
        // stage1query(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);
        stage1query2D(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints, device_queue);

        // lastError = cudaGetLastError();
        // if(lastError != cudaSuccess) printf("Error TRAVERSE: %s\n", cudaGetErrorString(lastError));
        syncError = cudaDeviceSynchronize();
        if(syncError != cudaSuccess) printf("Error TRAVERSE: %s\n", cudaGetErrorString(syncError));

        dtime = cast_t(tempo_t::now() - start).count();

        if(i%10 == 0)
            std::cout << " Partial " << i << " time elapsed: " << dtime << " ms\n";
        total += dtime;

    }

    uint32_t countMin=0;

#if DEVICE
    device_queue.memcpy(count_h, count, Ncells*sizeof(uint32_t)).wait();

    for(int i=0; i<Ncells; i++){
        if(count_h[i] != 0)
            countMin++;
    }

    free(count_h);
#else
    for(int i=0; i<Ncells; i++){
        if(count[i] != 0){
            countMin++;
            // printf("%u ", count[i]);
        }
    }
#endif

    std::cout << " Stage1 KERNEL time elapsed: " << total/n_tests << " ms\n";
    printf("Numero de minimos: %u\n", countMin);

    free(count, device_queue);

	return;
}
