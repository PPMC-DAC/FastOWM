#include <iostream>

#include <sycl_octree/source/octree_builder.h>

#include <basic/traverse_sycl.h>

void octree_traverse(std::string inputTXT, const uint32_t chunkDim, const uint32_t Wsize, const real_t Overlap)
{

#ifdef NVIDIA
    sycl::queue device_queue([](auto& d) 
        { return (d.get_platform().get_backend() == sycl::backend::ext_oneapi_cuda); });
#elif GPU
    sycl::queue device_queue(sycl::gpu_selector_v);
#elif CPU
    sycl::queue device_queue(sycl::cpu_selector_v);
#endif

    std::cout << "Device : " << device_queue.get_device().get_info<sycl::info::device::name>()  
    << " @ " << device_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << "Mhz (" << 
    device_queue.get_device().get_info<sycl::info::device::max_compute_units>() << " cores)" << std::endl;

	Octree_builder builder(inputTXT, chunkDim, device_queue);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.bintree.numObjects << " points; chunkDim: " << builder.bintree.leafSize << std::endl;
    std::cout << builder.bintree.numLeafs << " leaf nodes\n";
    std::cout << builder.bintree.numInternalNodes << " internal nodes\n";
    std::cout << builder.bintree.numNodes << " nodes\n\n";

    std::cout << "BBox: " << builder.bintree.BBox.upper.x << "," << builder.bintree.BBox.upper.y << " ";
    std::cout << builder.bintree.BBox.lower.x << "," << builder.bintree.BBox.lower.y << "\n\n";
  
    std::cout << "diffBox: " << builder.bintree.diffBox.x << "," << builder.bintree.diffBox.y << "\n\n";
#endif

    std::chrono::time_point<tempo_t> start = tempo_t::now();
    
    builder.build();

    std::chrono::time_point<tempo_t> end = tempo_t::now();
    // std::double_t dtime = cast_t(end - start).count();    
    //std::cout << "  CREATION takes: " << dtime << " ms\n";

    // uint32_t Wsize = 10;
    std::cout << "Wsize: " << Wsize << std::endl;
    // uint32_t Bsize = 20;
    // real_t Overlap = 0.99;
    std::cout << "Overlap: " << Overlap << std::endl;

    uint32_t Ncells;
    uint32_t nRows, nCols;
    real_t Width, High, Density;

    Width = builder.bintree.diffBox.x;
    High = builder.bintree.diffBox.y;
    // Average density of points
    Density = builder.bintree.numObjects/(Width*High);

    // Minimum number of points that should contain a Sliding Window so that its minimum is considered
    uint32_t minNumPoints = (uint32_t)(0.5*Density*Wsize*Wsize);

    real_t Displace = round2d(Wsize*(1-Overlap));

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
    int n_tests = 5;
    std::cout << "Performing " << n_tests << " tests\n";
#else
    int n_tests = 1;
#endif

    std::double_t total_s1 = 0.0, total_tree = 0.0;

    builder.reset();

    for(int i=0; i<n_tests; i++){

#if DEVICE
        device_queue.memset(count, 0u, Ncells*sizeof(uint32_t)).wait();
#else
        std::memset(count, 0u, Ncells*sizeof(uint32_t));
#endif

        start = tempo_t::now();

        builder.build();

        end = tempo_t::now();
        total_tree += cast_t(end - start).count();

        // stage1query(builder.node_list, builder.aabb_list, builder.ord_point_cloud, count,
        //         Wsize, Overlap, nCols, nRows, minNumPoints, builder.BBox, builder.diffBox, builder.numInternalNodes);
        // stage1query(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);

#ifdef CPU
        stage1query2DCPU(builder, count, Wsize, Overlap, nCols, 0u, nRows, minNumPoints);
#else
        auto e = stage1query2D(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints, device_queue);
        e.wait();
#endif
        // dtime = e.get_profiling_info<sycl::info::event_profiling::command_end>();
        total_s1 += cast_t(tempo_t::now() - end).count();

        builder.reset();
    }

    uint32_t countMin=0;

#if DEVICE
    start = tempo_t::now();
    device_queue.memcpy(count_h, count, Ncells*sizeof(uint32_t)).wait();

    // for(int i=0; i<Ncells; i++){
    //     if(count_h[i] != 0) countMin++;
    // }
    /* STAGE 2 */
    if(Overlap != 0.0){
        //qsort(count, Ncells, sizeof(uint32_t), &cmpfunc);
        std::sort(oneapi::dpl::execution::par_unseq,count_h, count_h+Ncells); //std::execution::par, std::execution::par_unseq,
        countMin = stage2CPU(Ncells, count_h);
        //printf("Numero de minimos STAGE2: %u\n", countMin);
    }
    double total_s2 = cast_t(tempo_t::now() - start).count();

    free(count_h);
#else
    start = tempo_t::now();
    // for(int i=0; i<Ncells; i++){
    //     if(count[i] != 0) countMin++;
    // }
    /* STAGE 2 */
    if(Overlap != 0.0){
        //qsort(count, Ncells, sizeof(uint32_t), &cmpfunc);
        std::sort(oneapi::dpl::execution::par_unseq,count, count+Ncells); //std::execution::par, std::execution::par_unseq,
        countMin = stage2CPU(Ncells, count);
        //printf("Numero de minimos STAGE2: %u\n", countMin);
    }
    double total_s2 = cast_t(tempo_t::now() - start).count();

#endif

    std::cout << " Tree Construction SYCL time elapased: " << total_tree/n_tests << " ms\n";
    std::cout << " Stage1 KERNEL SYCL time elapsed: " << total_s1/n_tests << " ms\n";
    std::cout << " Stage2 KERNEL SYCL time elapsed: " << total_s2 << " ms\n";
    std::cout << " Total KERNEL SYCL time elapsed: " << total_s1/n_tests + total_s2 << " ms\n";
    std::cout << " Total TIME (Tree+OWM) SYCL time elapsed: " << (total_tree + total_s1)/n_tests + total_s2<< " ms\n";
    printf("Numer of seed points: %u\n", countMin);

    free(count, device_queue);

	return;
}















void octree_traverse_heter(std::string inputTXT, const uint32_t chunkDim, const float factor, const uint32_t Wsize, const real_t Overlap)
{

// auto CUDASelector = [](sycl::device const &dev) {
//     if (dev.get_platform().get_backend() == sycl::backend::ext_oneapi_cuda) {
//       std::cout << " CUDA device found " << std::endl;
//       return 1;
//     } else {
//       return -1;
//     }
// };

// limit the maximum concurrency to that defined by the environment variable DPCPP_CPU_NUM_CUS for the TBB parallel_for in the CPU
// get the value of the environment variable
const char* num_cus = std::getenv("DPCPP_CPU_NUM_CUS");
int num_cpus = 0;
if(num_cus != NULL) {
    num_cpus = std::stoi(num_cus);
}else{
    num_cpus = 8;
}
tbb::global_control c(tbb::global_control::max_allowed_parallelism, num_cpus);

#ifdef NVIDIA
    sycl::queue device_queue([](auto& d) 
        { return (d.get_platform().get_backend() == sycl::backend::ext_oneapi_cuda); });
#elif GPU
    sycl::queue device_queue(sycl::gpu_selector_v);
#elif CPU
    sycl::queue device_queue(sycl::cpu_selector_v);
#endif

    std::cout << "Device : " << device_queue.get_device().get_info<sycl::info::device::name>()  
    << " @ " << device_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << "Mhz (" << 
    device_queue.get_device().get_info<sycl::info::device::max_compute_units>() << " cores)" << std::endl;

	Octree_builder builder(inputTXT, chunkDim, device_queue);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.bintree.numObjects << " points; chunkDim: " << builder.bintree.leafSize << std::endl;
    std::cout << builder.bintree.numLeafs << " leaf nodes\n";
    std::cout << builder.bintree.numInternalNodes << " internal nodes\n";
    std::cout << builder.bintree.numNodes << " nodes\n\n";

    std::cout << "BBox: " << builder.bintree.BBox.upper.x << "," << builder.bintree.BBox.upper.y << " ";
    std::cout << builder.bintree.BBox.lower.x << "," << builder.bintree.BBox.lower.y << "\n\n";
  
    std::cout << "diffBox: " << builder.bintree.diffBox.x << "," << builder.bintree.diffBox.y << "\n\n";
#endif

    std::chrono::time_point<tempo_t> start = tempo_t::now();
    
    builder.build();

    std::chrono::time_point<tempo_t> end = tempo_t::now();
    // std::double_t dtime = cast_t(tempo_t::now() - start).count();    
    // std::cout << "  CREATION takes: " << dtime << " ms\n";


    // uint32_t Wsize = 10;
    std::cout << "Wsize: " << Wsize << std::endl;
    // uint32_t Bsize = 20;
    // std::double_t Overlap = 0.99;
    std::cout << "Overlap: " << Overlap << std::endl;

    uint32_t Ncells;
    uint32_t nRows, nCols;
    std::double_t Width, High, Density;

    Width = builder.bintree.diffBox.x;
    High = builder.bintree.diffBox.y;
    // Densidad en puntos/m^2
    Density = builder.bintree.numObjects/(Width*High);

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
    // this vector will store the index of the selected minimums in the device memory
    uint32_t* count = (uint32_t*)sycl_builder::mallocWrap(Ncells*sizeof(uint32_t), device_queue);
    // the same vector but in the CPU
    uint32_t* count_cpu = (uint32_t*)std::malloc(Ncells*sizeof(uint32_t));

#if DEVICE
    // we need to copy the data from the device to the host
    uint32_t* count_h = (uint32_t*)std::malloc(Ncells*sizeof(uint32_t));
    octree_node* octree_h = (octree_node*)std::malloc(builder.m_node_count*sizeof(octree_node));
    aabb_t* aabb_h = (aabb_t*)std::malloc(builder.m_node_count*sizeof(aabb_t));
    point_t* points_h = (point_t*)std::malloc(builder.bintree.numObjects*sizeof(point_t));

    LBVHoct lbvh(
        octree_h,
        aabb_h,
        points_h
    );
#endif

#ifndef DEBUG
    int n_tests = 10;
    std::cout << "Performing " << n_tests << " tests (" << factor  << ")\n";
#else
    int n_tests = 1;
#endif

    std::double_t total_s1 = 0.0, total_tree = 0.0;
    // get the number of rows that the device will analyze
    uint32_t chunkRows = uint32_t(nRows*factor);

    builder.reset();

    for(int i=0; i<n_tests; i++){

#if DEVICE
        device_queue.memset(count, 0u, Ncells*sizeof(uint32_t)).wait();
        std::memset(count_cpu, 0u, Ncells*sizeof(uint32_t));
#else
        std::memset(count, 0u, Ncells*sizeof(uint32_t));
        std::memset(count_cpu, 0u, Ncells*sizeof(uint32_t));
#endif

        start = tempo_t::now();

        builder.build();

        end = tempo_t::now();
        total_tree += cast_t(end - start).count();

        // stage1query(builder.node_list, builder.aabb_list, builder.ord_point_cloud, count,
        //         Wsize, Overlap, nCols, nRows, minNumPoints, builder.BBox, builder.diffBox, builder.numInternalNodes);
        // stage1query(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);

        auto e = stage1query2D(builder, count, Wsize, Overlap, nCols, chunkRows, minNumPoints, device_queue);

#ifdef DEVICE
        // auto t_copy = tempo_t::now();

        device_queue.memcpy(octree_h, builder.m_octree, builder.m_node_count*sizeof(octree_node));
        device_queue.memcpy(aabb_h, builder.m_aabb, builder.m_node_count*sizeof(aabb_t));
        // we do not need a wait() here since the parallel_for, inside the stage1query2DCPU, will wait until the CPU finishes
        device_queue.memcpy(points_h, builder.bintree.ord_point_cloud, builder.bintree.numObjects*sizeof(point_t));

        // std::cout << "Copy time: " << cast_t(tempo_t::now() - t_copy).count() << "\n";

        stage1query2DCPU(lbvh, builder.bintree.BBox, count_cpu, Wsize, Overlap, nCols, std::make_pair(chunkRows, nRows), minNumPoints);

        // std::cout << "CPU time: " << cast_t(tempo_t::now() - t_copy).count() << "\n";

#else
        stage1query2DCPU(builder, count_cpu, Wsize, Overlap, nCols, std::make_pair(chunkRows, nRows), minNumPoints);
#endif

        e.wait();
        // dtime = e.get_profiling_info<sycl::info::event_profiling::command_end>();
        total_s1 += cast_t(tempo_t::now() - end).count();

        // device_queue.wait_and_throw();

        // if(i%10 == 0)
        //       std::cout << " Partial " << i << " time elapsed: " << dtime << " ms\n";
        // total += dtime;

        builder.reset();

    }

    uint32_t countMin=0;

#if DEVICE
    start = tempo_t::now();
    device_queue.memcpy(count_h, count, Ncells*sizeof(uint32_t)).wait();

    /* STAGE 2 */
    if(Overlap != 0.0){
        //qsort(count, Ncells, sizeof(uint32_t), &cmpfunc);
        std::sort(oneapi::dpl::execution::par_unseq, count_h, count_h+Ncells); //std::execution::par, std::execution::par_unseq,
        countMin = stage2CPU(Ncells, count_h);
        // printf("Numer of seed points: %u\n", countMin);
    }
    double total_s2 = cast_t(tempo_t::now() - start).count();

    free(count_h);
#else
    start = tempo_t::now();

    /* STAGE 2 */
    if(Overlap != 0.0){
        //qsort(count, Ncells, sizeof(uint32_t), &cmpfunc);
        std::sort(oneapi::dpl::execution::par_unseq, count, count+Ncells); //std::execution::par, std::execution::par_unseq,
        countMin = stage2CPU(Ncells, count);
        // printf("Numer of seed points: %u\n", countMin);
    }
    double total_s2 = cast_t(tempo_t::now() - start).count();

#endif

    std::cout << " Tree Construction SYCL time elapased: " << total_tree/n_tests << " ms\n";
    std::cout << " Stage1 KERNEL SYCL time elapsed: " << total_s1/n_tests << " ms\n";
    std::cout << " Stage2 KERNEL SYCL time elapsed: " << total_s2 << " ms\n";
    std::cout << " Total KERNEL SYCL time elapsed: " << total_s1/n_tests + total_s2 << " ms\n";
    std::cout << " Total TIME (Tree+OWM) SYCL time elapsed: " << (total_tree + total_s1)/n_tests + total_s2<< " ms\n";
    printf("Numer of seed points: %u\n", countMin);

    free(count, device_queue);
    free(count_cpu);

#ifdef DEVICE
    free(octree_h);
    free(aabb_h);
    free(points_h);
#endif

	return;
}
