#include <iostream>

#ifdef FEMU
#include <CL/sycl/INTEL/fpga_extensions.hpp>
#endif

#include <sycl_cuda_ordered/source/sycl_builder.h>

class CUDASelector : public sycl::device_selector {
public:
  int operator()(const sycl::device &Device) const override {
    using namespace sycl::info;

    const std::string DriverVersion = Device.get_info<device::driver_version>();

    if (Device.is_gpu() && (DriverVersion.find("CUDA") != std::string::npos)) {
      std::cout << "CUDA device found: ";
      return 1;
    };
    return -1;
  }
};


void sycl_test(std::string inputTXT, uint32_t chunkDim)
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

    double time = 0.0;
#ifndef DEBUG
    int n_tests = 100;
    std::cout << "Performing " << n_tests+1 << " tests:\n";
#else
    int n_tests = 1;
    std::cout << "Performing " << n_tests+1 << " test:\n";
#endif

    for (int i = 0; i <= n_tests; i++)
    {
        std::chrono::time_point<tempo_t> start = tempo_t::now();
      
        builder.build();

        // cudaError_t lastError = cudaGetLastError();
        // if(lastError != cudaSuccess) printf("Error RESET: %s\n", cudaGetErrorString(lastError));
        cudaError_t syncError = cudaDeviceSynchronize();
        if(syncError != cudaSuccess) printf("Sync error RESET: %s\n", cudaGetErrorString(syncError));
    
        double dtime = cast_t(tempo_t::now() - start).count();
        
        if(i)
          time += dtime;

#ifdef DEBUG
        std::cout << "Test: " << dtime << " ms\n";
#endif

        builder.reset();
    }

	time /= (float)(n_tests);
    
    std::cout << "  CREATION takes: " << time << " ms\n";

	return;
}
