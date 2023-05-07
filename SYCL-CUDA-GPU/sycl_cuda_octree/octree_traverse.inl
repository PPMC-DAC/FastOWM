#include <iostream>

#include <sycl_cuda_octree/cuda/octree_builder.h>

#include <basic/traverse_sycl_cuda_octree.h>









void octree_traverse(std::string inputTXT, const uint32_t chunkDim)
{

    cl::sycl::queue device_queue(CUDASelector().select_device());

    std::cout << device_queue.get_device().get_info<sycl::info::device::name>()  
    << " @ " << device_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << "Mhz (" << 
    device_queue.get_device().get_info<sycl::info::device::max_compute_units>() << " cores)" << std::endl;

    uint32_t numObjects = 0;
    whole_t BBox;

#ifndef CHECK
    readHeader(inputTXT, BBox, numObjects);
#else
    inputTXT = "CHECK";
    chunkDim = 2;
    numObjects = 32;
#endif

    uint32_t numLeafs = (uint32_t)((numObjects-1)/chunkDim) + 1u;
    uint32_t numInternalNodes = numLeafs - 1u;
    uint32_t numNodes = 2*numLeafs - 1u;

    // Numero minimo de nodos en bloques de 8 que contienen a todas las hojas
    uint32_t aux = uint32_t((numLeafs-1)/8) + 1u;
    // pero necesito todos los nodos intermedios también
    int i=0;
    uint32_t maxNumTasks = 0;
    uint32_t pot = 0;
    for(; pot<aux; i++) {
        pot = uint32_t(pow(8,i));
        maxNumTasks += pot;
    }

    /* Procesar un nodo hoja no va a ser una tarea para nosotros, 
    pero sí que va a ser un nuevo nodo del octree*/

    uint32_t maxNumNodes = maxNumTasks + numLeafs;

#ifdef DEBUG
    printf("Necesito %d niveles en el arbol\n", i);
    printf("Voy a tener como máximo %u nodos\n", maxNumNodes);
#endif

    // SYCL_builder builder(inputTXT, chunkDim, device_queue, numObjects, numNodes, numInternalNodes);

	Octree_builder builder(inputTXT, chunkDim, device_queue, numObjects, numNodes, numInternalNodes, maxNumTasks, maxNumNodes);

    std::chrono::time_point<tempo_t> start = tempo_t::now();
    
    builder.build();

    cudaError_t syncError = cudaDeviceSynchronize();
    if(syncError != cudaSuccess) printf("Sync error TRAVERSE: %s\n", cudaGetErrorString(syncError));

    double dtime = cast_t(tempo_t::now() - start).count();    
    std::cout << "  CREATION takes: " << dtime << " ms\n";


    uint32_t Wsize = 10;
    // uint32_t Bsize = 20;
    double Overlap = 0.8;

    uint32_t Ncells;
    uint32_t nRows, nCols;
    double Width, High, Density;

    uint32_t countMin;

    Width = builder.bintree.diffBox.x;
    High = builder.bintree.diffBox.y;
    // Densidad en puntos/m^2
    Density = builder.bintree.numObjects/(Width*High);

    // El numero minimo sera la mitad del numero de puntos medio por celda
    uint32_t minNumPoints = (uint32_t)(0.5*Density*Wsize*Wsize);

    double Displace = round2d(Wsize*(1-Overlap));

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
    printf("%d,%d\n", nCols, nRows);
#endif

    /*Este es donde ya lanzo un kernel CUDA de búsqueda*/
    // uint32_t* count = NULL;

    // cudaMallocManaged((void**)&count, Ncells*sizeof(uint32_t));

    // sycl::host_accessor accN(builder.b_octree, sycl::read_only);
    // sycl::host_accessor accA(builder.b_aabb, sycl::read_only);
    // sycl::host_accessor accP(builder.bintree.ord_point_cloud, sycl::read_only);

    // LBVHoctCPU lbvh(
    //     accN,
    //     accA,
    //     accP
    // );

    sycl::buffer<uint32_t, 1> b_count(Ncells);

    // uint32_t* cpu_count = (uint32_t*)std::malloc(Ncells * sizeof(uint32_t));

    // std::memset(cpu_count, 0u, Ncells * sizeof(uint32_t));

#ifndef DEBUG
    int n_tests = 50;
    std::cout << "Performing " << n_tests+1 << " tests\n";
#else
    int n_tests = 1;
    std::cout << "Performing " << n_tests+1 << " tests:\n";
#endif
    double total = 0.0;

    builder.reset();

    for(int i=0; i<=n_tests; i++){
        // cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);

        start = tempo_t::now();

        builder.build();

        // stage1query(builder, b_count, Wsize, Overlap, nCols, nRows, minNumPoints, device_queue);
        stage1query2D(builder, b_count, Wsize, Overlap, nCols, nRows, minNumPoints, device_queue);
        // stage1query2DCPU(lbvh, builder.bintree.BBox, cpu_count, Wsize, Overlap, 0u, nCols, nRows, minNumPoints);

        cudaError_t syncError = cudaDeviceSynchronize();
        if(syncError != cudaSuccess) printf("Sync error TRAVERSE: %s\n", cudaGetErrorString(syncError));
        // device_queue.wait();

        dtime = cast_t(tempo_t::now() - start).count();
        if(i%10 == 0)
            std::cout << " Partial " << i << " time elapsed: " << dtime << " ms\n";

        builder.reset();

        if(i)
            total += dtime;

        // cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);
    }

    {
        sycl::host_accessor h_count(b_count);
        countMin=0;
        for(int i=0; i<Ncells; i++){
            if(h_count[i] != 0)
                countMin++;
            h_count[i] = 0u;
        }
        // for(int i=0; i<Ncells; i++){
        //     if(cpu_count[i] != 0)
        //         countMin++;
        //     cpu_count[i] = 0u;
        // }
    }    

    std::cout << " Stage1 KERNEL CUDA time elapsed: " << total/n_tests << " ms\n";
    printf("Numero de minimos: %u\n", countMin);

    // cudaFree(count);
    // free(cpu_count);

	return;
}












void octree_traverse_heter(std::string inputTXT, const uint32_t chunkDim, const float factor)
{

    cl::sycl::queue device_queue(CUDASelector().select_device());

    std::cout << device_queue.get_device().get_info<sycl::info::device::name>()  
    << " @ " << device_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << "Mhz (" << 
    device_queue.get_device().get_info<sycl::info::device::max_compute_units>() << " cores)" << std::endl;

    uint32_t numObjects = 0;
    whole_t BBox;

#ifndef CHECK
    readHeader(inputTXT, BBox, numObjects);
#else
    inputTXT = "CHECK";
    chunkDim = 2;
    numObjects = 32;
#endif

    uint32_t numLeafs = (uint32_t)((numObjects-1)/chunkDim) + 1u;
    uint32_t numInternalNodes = numLeafs - 1u;
    uint32_t numNodes = 2*numLeafs - 1u;

    // Numero minimo de nodos en bloques de 8 que contienen a todas las hojas
    uint32_t aux = uint32_t((numLeafs-1)/8) + 1u;
    // pero necesito todos los nodos intermedios también
    int i=0;
    uint32_t maxNumTasks = 0;
    uint32_t pot = 0;
    for(; pot<aux; i++) {
        pot = uint32_t(pow(8,i));
        maxNumTasks += pot;
    }

    /* Procesar un nodo hoja no va a ser una tarea para nosotros, 
    pero sí que va a ser un nuevo nodo del octree*/

    uint32_t maxNumNodes = maxNumTasks + numLeafs;

#ifdef DEBUG
    printf("Necesito %d niveles en el arbol\n", i);
    printf("Voy a tener como máximo %u nodos\n", maxNumNodes);
#endif

    // SYCL_builder builder(inputTXT, chunkDim, device_queue, numObjects, numNodes, numInternalNodes);

	Octree_builder builder(inputTXT, chunkDim, device_queue, numObjects, numNodes, numInternalNodes, maxNumTasks, maxNumNodes);

    std::chrono::time_point<tempo_t> start = tempo_t::now();
    
    builder.build();

    cudaError_t syncError = cudaDeviceSynchronize();
    if(syncError != cudaSuccess) printf("Sync error TRAVERSE: %s\n", cudaGetErrorString(syncError));

    double dtime = cast_t(tempo_t::now() - start).count();    
    std::cout << "  CREATION takes: " << dtime << " ms\n";


    uint32_t Wsize = 10;
    // uint32_t Bsize = 20;
    double Overlap = 0.8;

    uint32_t Ncells;
    uint32_t nRows, nCols;
    double Width, High, Density;

    uint32_t countMin;

    Width = builder.bintree.diffBox.x;
    High = builder.bintree.diffBox.y;
    // Densidad en puntos/m^2
    Density = builder.bintree.numObjects/(Width*High);

    // El numero minimo sera la mitad del numero de puntos medio por celda
    uint32_t minNumPoints = (uint32_t)(0.5*Density*Wsize*Wsize);

    double Displace = round2d(Wsize*(1-Overlap));

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
    printf("%d,%d\n", nCols, nRows);
    printf("offload factor: %f\n", factor);
#endif

    sycl::host_accessor accN(builder.b_octree, sycl::read_only);
    sycl::host_accessor accA(builder.b_aabb, sycl::read_only);
    sycl::host_accessor accP(builder.bintree.ord_point_cloud, sycl::read_only);

    LBVHoctCPU lbvh(
        accN,
        accA,
        accP
    );

    sycl::buffer<uint32_t, 1> b_count(Ncells);

    {
        sycl::host_accessor h_count(b_count);
        for(int i=0; i<Ncells; i++){
            h_count[i] = 0u;
        }
    }

    uint32_t* cpu_count = (uint32_t*)std::malloc(Ncells * sizeof(uint32_t));

    std::memset(cpu_count, 0u, Ncells * sizeof(uint32_t));

#ifndef DEBUG
    int n_tests = 50;
    std::cout << "Performing " << n_tests+1 << " tests (" << factor  << ")\n";
#else
    int n_tests = 1;
    std::cout << "Performing " << n_tests+1 << " tests:\n";
#endif
    double total = 0.0;

    uint32_t wCols = uint32_t(nCols*factor);

    // builder.reset();

    for(int i=0; i<=n_tests; i++){
        // cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);

        start = tempo_t::now();

        // builder.build();

        stage1query2D(builder, b_count, Wsize, Overlap, wCols, nRows, minNumPoints, device_queue);
        // stage1query2DCPU(builder, cpu_count, Wsize, Overlap, wCols, nCols, nRows, minNumPoints);
        stage1query2DCPU(lbvh, builder.bintree.BBox, cpu_count, Wsize, Overlap, wCols, nCols, nRows, minNumPoints);

        // {
        //     sycl::host_accessor accN(builder.b_octree, sycl::read_only);
        //     sycl::host_accessor accA(builder.b_aabb, sycl::read_only);
        //     sycl::host_accessor accP(builder.bintree.ord_point_cloud, sycl::read_only);

        //     LBVHoctCPU lbvh(
        //         accN,
        //         accA,
        //         accP
        //     );

        //     stage1query2DCPU(lbvh, builder.bintree.BBox, cpu_count, Wsize, Overlap, wCols, nCols, nRows, minNumPoints);
        // }

        cudaError_t syncError = cudaDeviceSynchronize();
        if(syncError != cudaSuccess) printf("Sync error TRAVERSE: %s\n", cudaGetErrorString(syncError));
        // device_queue.wait();

        dtime = cast_t(tempo_t::now() - start).count();
        if(i%10 == 0)
            std::cout << " Partial " << i << " time elapsed: " << dtime << " ms\n";

        // builder.reset();

        if(i)
            total += dtime;

        // cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);
    }

    {
        sycl::host_accessor h_count(b_count);
        countMin=0;
        for(int i=0; i<Ncells; i++){
            if(h_count[i] != 0)
                countMin++;
            h_count[i] = 0u;
        }
        for(int i=0; i<Ncells; i++){
            if(cpu_count[i] != 0)
                countMin++;
            cpu_count[i] = 0u;
        }
    }    

    std::cout << " Stage1 KERNEL CUDA time elapsed: " << total/n_tests << " ms\n";
    printf("Numero de minimos: %u\n", countMin);

    // cudaFree(count);
    free(cpu_count);

	return;
}
