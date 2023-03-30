#include <iostream>

#include <octree/cuda/octree_builder.h>

#include <basic/traverse.h>

void octree_traverse(std::string inputTXT, const uint32_t chunkDim)
{

	Octree_builder builder(inputTXT, chunkDim);

    std::chrono::time_point<tempo_t> start = tempo_t::now();
    builder.build();
    std::chrono::time_point<tempo_t> end = tempo_t::now();

    double dtime = cast_t(end - start).count();    
   
    uint32_t Wsize = 10;
    // uint32_t Bsize = 20;
    double Overlap = 0.8;

    uint32_t Ncells;
    uint32_t nRows, nCols;
    double Width, High, Density;

    uint32_t countMin;

    Width = builder.bintree.diffBox.x;
    High = builder.bintree.diffBox.y;
    // Average density of points
    Density = builder.bintree.numObjects/(Width*High);

    // Minimum number of points that should contain a Sliding Window so that its minimum is considered
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

    uint32_t* count = NULL;

    cudaMallocManaged((void**)&count, Ncells*sizeof(uint32_t));

#ifndef DEBUG
    int n_tests = 5;
#else
    int n_tests = 1;
#endif
    double total_s1{0.0}, total_s2{0.0}, total_tree{0.0};

    builder.reset();

    for(int i=0; i<n_tests; i++){
        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);

        start = tempo_t::now();

        builder.build();
        end = tempo_t::now();
        total_tree += cast_t(end - start).count();

        // stage1query(builder.node_list, builder.aabb_list, builder.ord_point_cloud, count,
        //         Wsize, Overlap, nCols, nRows, minNumPoints, builder.BBox, builder.diffBox, builder.numInternalNodes);
        // stage1query(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);
        stage1query2D(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);

        cudaDeviceSynchronize();

        total_s1 += cast_t(tempo_t::now() - end).count();

        start = tempo_t::now();
        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);

        for(int i=0; i<Ncells; i++){
            if(count[i] != 0xFFFFFFFFu)
                countMin++;
        }
        //printf("Numero de minimos STAGE1: %u\n", countMin);

        /* STAGE 2 */

        if(Overlap != 0.0){
            //qsort(count, Ncells, sizeof(uint32_t), &cmpfunc);
            std::sort(count, count+Ncells); //std::execution::par, std::execution::par_unseq,
            countMin = stage2CPU(Ncells, count);
            //printf("Numero de minimos STAGE2: %u\n", countMin);
        }
        total_s2 += cast_t(tempo_t::now() - start).count();

        // Fichero de salida
        // if(save_file("prueba_INAER_octree.xyz", count, countMin, builder.bintree.ord_point_cloud) < 0){
        //     printf("Unable to create results file!\n");
        // }
    }
    std::cout << " Tree Construction CUDA time elapased: " << total_tree/n_tests << " ms\n";
    std::cout << " Stage1 KERNEL CUDA time elapsed: " << total_s1/n_tests << " ms\n";
    std::cout << " Stage2 KERNEL CUDA time elapsed: " << total_s2/n_tests << " ms\n";
    std::cout << " Total KERNEL CUDA time elapsed: " << (total_s1+total_s2)/n_tests << " ms\n";
    std::cout << " Total TIME (Tree+OWM) CUDA time elapsed: " << total_tree + (total_s1+total_s2)/n_tests << " ms\n";
    printf("Numer of seed points: %u\n", countMin);

    cudaFree(count);

	return;
}


void octree_traverse_heter(std::string inputTXT, const uint32_t chunkDim, const float factor)
{

	Octree_builder builder(inputTXT, chunkDim);

    std::chrono::time_point<tempo_t> start = tempo_t::now();
    builder.build();
    std::chrono::time_point<tempo_t> end = tempo_t::now();
    double dtime = cast_t(end - start).count();    
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

    uint32_t wCols = uint32_t(nCols*factor);

    uint32_t* count = NULL;
    uint32_t* cpu_count = NULL;

    cudaMallocManaged((void**)&count, Ncells*sizeof(uint32_t));
    cudaMallocHost((void**)&cpu_count, Ncells*sizeof(uint32_t));

    // cudaMemPrefetchAsync(builder.b_octree, builder.maxNumNodes*sizeof(octree_node), cudaCpuDeviceId);
    // cudaMemPrefetchAsync(builder.b_aabb, builder.maxNumNodes*sizeof(aabb_t), cudaCpuDeviceId);
    // cudaMemPrefetchAsync(builder.bintree.ord_point_cloud, builder.bintree.numObjects*sizeof(point_t), cudaCpuDeviceId);

    // cudaMemset(count, 0u, Ncells*sizeof(uint32_t));
    // std::memset(cpu_count, 0u, Ncells*sizeof(uint32_t));

#ifndef DEBUG
    int n_tests = 5;
    std::cout << "Performing " << n_tests << " tests (" << chunkDim << ", " << factor  << ")\n";
#else
    int n_tests = 1;
#endif
    double total_s1{0.0}, total_tree{0.0};

    builder.reset();

    for(int i=0; i<=n_tests; i++){

        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);

        start = tempo_t::now();
        builder.build();
        end = tempo_t::now();
        total_tree += cast_t(end - start).count();

        // stage1query(builder.node_list, builder.aabb_list, builder.ord_point_cloud, count,
        //         Wsize, Overlap, nCols, nRows, minNumPoints, builder.BBox, builder.diffBox, builder.numInternalNodes);
        // stage1query(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);
        stage1query2D(builder, count, Wsize, Overlap, wCols, nRows, minNumPoints);

        stage1query2DCPU(builder, cpu_count, Wsize, Overlap, wCols, nCols, nRows, minNumPoints);

        cudaDeviceSynchronize();

        total_s1 += cast_t(tempo_t::now() - end).count();

        builder.reset();

        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);

        countMin=0;
        for(int i=0; i<Ncells; i++){
            if(count[i] != 0)
                countMin++;
            count[i] = 0u;
        }
        for(int i=0; i<Ncells; i++){
            if(cpu_count[i] != 0)
                countMin++;
            cpu_count[i] = 0u;
        }
    }
    std::cout << " Tree Construction CUDA time elapased: " << total_tree/n_tests << " ms\n";
    std::cout << " Stage1 KERNEL CUDA time elapsed: " << total_s1/n_tests << " ms\n";

    printf("Number of minima: %u\n", countMin);

    cudaFree(count);
    cudaFreeHost(cpu_count);

	return;
}
