#include <iostream>

#ifdef FEMU
#include <CL/sycl/INTEL/fpga_extensions.hpp>
#endif

#include <sycl_ordered/source/sycl_builder.h>

void sycl_test(std::string inputTXT, uint32_t chunkDim)
{

#ifdef CHECK
    inputTXT = "CHECK";
    chunkDim = 2;
#endif

#ifdef FEMU
    cl::sycl::INTEL::fpga_emulator_selector device_selector{};
#elif GPU
    cl::sycl::gpu_selector device_selector{};
#elif CPU
    cl::sycl::cpu_selector device_selector{};
#endif

    cl::sycl::queue device_queue(device_selector);

    std::cout << "Device : " << device_queue.get_device().get_info<sycl::info::device::name>()  
    << " @ " << device_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << "Mhz (" << 
    device_queue.get_device().get_info<sycl::info::device::max_compute_units>() << " cores)" << std::endl;

    SYCL_builder builder(inputTXT, chunkDim, device_queue);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.numObjects << " points; chunkDim: " << builder.leafSize << std::endl;
    std::cout << builder.numLeafs << " leaf nodes\n";
    std::cout << builder.numInternalNodes << " internal nodes\n";
    std::cout << builder.numNodes << " nodes\n\n";

    std::cout << "BBox: " << builder.BBox.upper.x << "," << builder.BBox.upper.y << " ";
    std::cout << builder.BBox.lower.x << "," << builder.BBox.lower.y << "\n\n";
  
    std::cout << "diffBox: " << builder.diffBox.x << "," << builder.diffBox.y << "\n\n";
#endif

    std::double_t time = 0.0;
#ifndef DEBUG
    int n_tests = 100;
    std::cout << "Performing " << n_tests << " tests:\n";
#else
    int n_tests = 2;
    std::cout << "Performing " << n_tests << " test:\n";
#endif

    for (int i = 0; i < n_tests; i++)
    {
        std::chrono::time_point<tempo_t> start = tempo_t::now();
      
        builder.build();
    
        std::double_t dtime = cast_t(tempo_t::now() - start).count();
        
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
