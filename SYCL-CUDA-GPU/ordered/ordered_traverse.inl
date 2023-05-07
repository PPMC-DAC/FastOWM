#include <iostream>

#include <ordered/cuda/ordered_builder.h>

#include <basic/traverse.h>

#include <basic/qtree.h>


void ordered_traverse(std::string inputTXT, const uint32_t chunkDim)
{

	Ordered_builder builder(inputTXT, chunkDim);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.numObjects << " points; chunkDim: " << builder.leafSize;
    std::cout << "; threadsPerBlock: " << builder.threadsPerBlock << "; numBlocks: " << builder.numberOfBlocks << std::endl;
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
    lastError = cudaDeviceSynchronize();
    if(lastError != cudaSuccess) printf("Error BUILD SYNC: %s\n", cudaGetErrorString(lastError));

    double dtime = cast_t(tempo_t::now() - start).count();    
    std::cout << "  CREATION takes: " << dtime << " ms\n";

    double build_time = dtime;

    uint32_t Wsize = 10;
    uint32_t Bsize = 20;
    double Overlap = 0.8;

    uint32_t Ncells, Ngrids;
    uint32_t nRows, nCols, nRowsg, nColsg;
    double Width, High, Density;

    uint32_t countMin;

    Width = builder.diffBox.x;
    High = builder.diffBox.y;
    // Densidad en puntos/m^2
    Density = builder.numObjects/(Width*High);

    // El numero minimo sera la mitad del numero de puntos medio por celda
    uint32_t minNumPoints = (uint32_t)(0.5*Density*Wsize*Wsize);

    double Displace = round2d(Wsize*(1.-Overlap));

    // Stage 1 parameters
    if(Overlap > 0.0) {
        nCols=(uint32_t)(round((Width+2*Wsize*Overlap)/Displace))-1;
        nRows=(uint32_t)(round((High+2*Wsize*Overlap)/Displace))-1;
    } else {
        nCols=(uint32_t)floor(Width/Wsize)+1;
        nRows=(uint32_t)floor(High/Wsize)+1;
    }
    Ncells = nCols*nRows;
#ifdef DEBUG
    printf("%d,%d\n", nCols, nRows);
#endif

    // Stage 3 parameters
    nColsg=(uint32_t)floor(Width/Bsize)+1;
    nRowsg=(uint32_t)floor(High/Bsize)+1;
    Ngrids = nColsg*nRowsg;
#ifdef DEBUG
    printf("%d,%d\n", nColsg, nRowsg);
#endif

    /*Reservo el espacio para el árbol de S3, que va a tener como máximo "Ncells" puntos*/
    // Ordered_builder builder_grid(builder.BBox, Ncells, chunkDim);

    /*Este es donde ya lanzo un kernel CUDA de búsqueda*/
    uint32_t* count = NULL;
    uint32_t* countStage3 = NULL;

    cudaMallocManaged((void**)&count, Ncells*sizeof(uint32_t));
    cudaMallocManaged((void**)&countStage3, Ngrids*sizeof(uint32_t));

#ifndef DEBUG
    int n_tests = 100;
#else
    int n_tests = 1;
#endif
    double total = 0.0;

    for(int i=0; i<n_tests; i++){

        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);

        start = tempo_t::now();

        /* STAGE 1 */

        // stage1query(builder.node_list, builder.aabb_list, builder.ord_point_cloud, count,
        //         Wsize, Overlap, nCols, nRows, minNumPoints, builder.BBox, builder.diffBox, builder.numInternalNodes);
        // stage1query(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);
        stage1query2D(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);

        dtime = cast_t(tempo_t::now() - start).count();
        if(i%10 == 0)
            std::cout << " Partial S1 " << i << " time elapsed: " << dtime << " ms\n";
        total += dtime;

        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);

#ifdef DEBUG

        for(int j=0; j<Ncells; j++){
            if(count[j] != 0xFFFFFFFFu)
                countMin++;
        }
        printf("Numero de minimos STAGE1: %u\n", countMin);

#endif

        /* STAGE 2 */

        start = tempo_t::now();

        if(Overlap != 0.0){
            qsort(count, Ncells, sizeof(uint32_t), &cmpfunc);
            countMin = stage2CPU(Ncells, count);
        }

        dtime = cast_t(tempo_t::now() - start).count();
        if(i%10 == 0)
            std::cout << " Partial S2 " << i << " time elapsed: " << dtime << " ms\n";
        total += dtime;

#ifdef DEBUG
        printf("Numero de minimos STAGE2: %u\n", countMin);
#endif

        /* STAGE 3 */

        start = tempo_t::now();

        // Aplico la malla para ver si hay zonas sin puntos
        if(Bsize > 0){ 

            cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);
            cudaMemPrefetchAsync(countStage3, Ngrids*sizeof(uint32_t), 0);

            Ordered_builder builder_grid( builder.ord_point_cloud, builder.BBox, count, countMin, chunkDim );

            builder_grid.build();

            lastError = cudaGetLastError();
            if(lastError != cudaSuccess) printf("Error BUILD GRID: %s\n", cudaGetErrorString(lastError));
            lastError = cudaDeviceSynchronize();
            if(lastError != cudaSuccess) printf("Error BUILD GRID SYNC: %s\n", cudaGetErrorString(lastError));

            stage3query2D(builder, builder_grid, countStage3, Bsize, nColsg, nRowsg, minNumPoints);

            builder_grid.reset();

        }// Bsize

        dtime = cast_t(tempo_t::now() - start).count();
        if(i%10 == 0)
            std::cout << " Partial S3 " << i << " time elapsed: " << dtime << " ms\n";
        total += dtime;

        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);
        cudaMemPrefetchAsync(countStage3, Ngrids*sizeof(uint32_t), cudaCpuDeviceId);

        if(i != n_tests-1){
            countMin=0;
            for(int j=0; j<Ncells; j++){
                count[j] = 0xFFFFFFFFu;
            }
            // countMinStage3=0;
            for(int j=0; j<Ngrids; j++){
                countStage3[j] = 0xFFFFFFFFu;
            }
        }

    }
    std::cout << " CUDA SEARCH time elapsed: " << total/n_tests << " ms\n";

    std::vector<uint32_t> ids;

    for(int j=0; j<countMin; j++){
        if(count[j] != 0xFFFFFFFFu)
            ids.push_back(count[j]);
    }

    uint32_t countMinStage3 = 0u;
    for(int j=0; j<Ngrids; j++){
        if(countStage3[j] != 0xFFFFFFFFu){
            ids.push_back(countStage3[j]);
            countMinStage3++;
        }
    }
    printf("Numero de minimos STAGE3: %u\n", countMinStage3);
    printf("Número total de mínimos: %zu\n", ids.size());

    // Fichero de salida
    if(save_file("prueba_INAER_ordered.xyz", count, countMin, builder.ord_point_cloud) < 0){
        printf("Unable to create results file!\n");
    }
    if(save_file_append("prueba_INAER_ordered.xyz", countStage3, Ngrids, builder.ord_point_cloud) < 0){
        printf("Unable to create results file!\n");
    }

    std::string gold_results = obtain_gold(inputTXT);

    double checked_rate = check_results(gold_results, ids, builder.ord_point_cloud, Displace);
    if (checked_rate < 0) {
        printf("Unable to check results\n");
    }

    if(save_time("resultados_ORDERED.csv", inputTXT, chunkDim, builder.threadsPerBlock, builder.numberOfBlocks, build_time, total/n_tests, checked_rate) < 0){
        printf("Unable to create results file!\n");
    }

    builder.reset();
    
    cudaFree(count);
    cudaFree(countStage3);

	return;
}









void ordered_traverse_S3CPU(std::string inputTXT, const uint32_t chunkDim)
{

	Ordered_builder builder(inputTXT, chunkDim);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.numObjects << " points; chunkDim: " << builder.leafSize;
    std::cout << "; threadsPerBlock: " << builder.threadsPerBlock << "; numBlocks: " << builder.numberOfBlocks << std::endl;
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
    cudaDeviceSynchronize();

    double dtime = cast_t(tempo_t::now() - start).count();    
    std::cout << "  CREATION takes: " << dtime << " ms\n";


    uint32_t Wsize = 12;
    uint32_t Bsize = 20;
    double Overlap = 0.5;

    uint32_t Ncells, Ngrids;
    uint32_t nRows, nCols, nRowsg, nColsg;
    double Width, High, Density;

    uint32_t countMin;

    Width = builder.diffBox.x;
    High = builder.diffBox.y;
    // Densidad en puntos/m^2
    Density = builder.numObjects/(Width*High);

    // El numero minimo sera la mitad del numero de puntos medio por celda
    uint32_t minNumPoints = (uint32_t)(0.5*Density*Wsize*Wsize);

    double Displace = round2d(Wsize*(1-Overlap));

    // Stage 1 parameters
    if(Overlap > 0.0) {
        nCols=(uint32_t)(round((Width+2*Wsize*Overlap)/Displace))-1;
        nRows=(uint32_t)(round((High+2*Wsize*Overlap)/Displace))-1;
    } else {
        nCols=(uint32_t)floor(Width/Wsize)+1;
        nRows=(uint32_t)floor(High/Wsize)+1;
    }
    Ncells = nCols*nRows;
#ifdef DEBUG
    printf("%d,%d\n", nCols, nRows);
#endif

    // Stage 3 parameters
    nColsg=(uint32_t)floor(Width/Bsize)+1;
    nRowsg=(uint32_t)floor(High/Bsize)+1;
    Ngrids = nColsg*nRowsg;
#ifdef DEBUG
    printf("%d,%d\n", nColsg, nRowsg);
#endif

    real_v center, radius, min, max;
    float maxRadius = 0.0;
    min = {builder.BBox.lower.x, builder.BBox.lower.y};
    max = {builder.BBox.upper.x, builder.BBox.upper.y};
    radius = getRadius(min, max, &maxRadius);
    center = getCenter(min, radius);

    /*Este es donde ya lanzo un kernel CUDA de búsqueda*/
    uint32_t* count = NULL;
    uint32_t* countStage3 = NULL;

    cudaMallocManaged((void**)&count, Ncells*sizeof(uint32_t));
    cudaMallocManaged((void**)&countStage3, Ngrids*sizeof(uint32_t));

#ifndef DEBUG
    int n_tests = 40;
#else
    int n_tests = 1;
#endif
    double total = 0.0;

    for(int i=0; i<n_tests; i++){
        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);

        start = tempo_t::now();

        /* STAGE 1 */

        // stage1query(builder.node_list, builder.aabb_list, builder.ord_point_cloud, count,
        //         Wsize, Overlap, nCols, nRows, minNumPoints, builder.BBox, builder.diffBox, builder.numInternalNodes);
        // stage1query(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);
        stage1query2D(builder, count, Wsize, Overlap, nCols, nRows, minNumPoints);

        dtime = cast_t(tempo_t::now() - start).count();
        if(i%10 == 0)
            std::cout << " Partial S1 " << i << " time elapsed: " << dtime << " ms\n";
        total += dtime;

        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);

#ifdef DEBUG

        for(int j=0; j<Ncells; j++){
            if(count[j] != 0xFFFFFFFFu)
                countMin++;
        }
        printf("Numero de minimos STAGE1: %u\n", countMin);

#endif

        /* STAGE 2 */

        start = tempo_t::now();

        if(Overlap != 0.0){
            qsort(count, Ncells, sizeof(uint32_t), &cmpfunc);
            countMin = stage2CPU(Ncells, count);
        }

        dtime = cast_t(tempo_t::now() - start).count();
        if(i%10 == 0)
            std::cout << " Partial S2 " << i << " time elapsed: " << dtime << " ms\n";
        total += dtime;

#ifdef DEBUG
        printf("Numero de minimos STAGE2: %u\n", countMin);
#endif

        /* STAGE 3 */

        start = tempo_t::now();

        // cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), 0);
        // cudaMemPrefetchAsync(countStage3, Ngrids*sizeof(uint32_t), 0);

        // // Aplico la malla para ver si hay zonas sin puntos
        // if(Bsize > 0){ 

        //     Ordered_builder builder_grid(builder.ord_point_cloud, builder.numObjects, builder.BBox, count, countMin, chunkDim);

        //     builder_grid.build();

        //     stage3query2D(builder, builder_grid, countStage3, Bsize, nColsg, nRowsg, minNumPoints);

        //     builder_grid.reset();

        // }// Bsize

        Qtree grid = createQtree(NULL, center, maxRadius);

        if(Overlap > 0.0){
            for(uint32_t j = 0; j < countMin; j++){
                // point_t p = builder.ord_point_cloud[count[j]];
                // Lpoint lp = {j,p.x,p.y,p.z};
                insertPoint2(&builder.ord_point_cloud[count[j]], grid, 128);
            }
        } else { // because data is all out of order; NO stage2
            // for(int j = 0; j < Ncells; j++)
            //     if(count[j] != 0xFFFFFFFFu) insertPoint2(&builder.ord_point_cloud[count[j]], grid, 128);
        }

        /*AUN FALTA LA OBTENECION DE LOS PUNTOS, PERO YA TARDA MÁS, ASÍ QUE LO VAMOS A DEJAR AQUÍ*/

        deleteQtree(grid);

        delete(grid);

        dtime = cast_t(tempo_t::now() - start).count();
        if(i%10 == 0)
            std::cout << " Partial S3 " << i << " time elapsed: " << dtime << " ms\n";
        total += dtime;

#ifdef DEBUG
        uint32_t countMinStage3 = 0u;
        for(int j=0; j<Ngrids; j++){
            if(countStage3[j] != 0xFFFFFFFFu)
                countMinStage3++;
        }
        printf("Numero de minimos STAGE3: %u\n", countMinStage3);

        // Fichero de salida
        if(save_file("prueba_INAER_ordered.xyz", count, countMin, builder.ord_point_cloud) < 0){
            printf("Unable to create results file!\n");
        }
        if(save_file_append("prueba_INAER_ordered.xyz", countStage3, Ngrids, builder.ord_point_cloud) < 0){
            printf("Unable to create results file!\n");
        }
#endif

        // builder.reset();

        cudaMemPrefetchAsync(count, Ncells*sizeof(uint32_t), cudaCpuDeviceId);
        cudaMemPrefetchAsync(countStage3, Ngrids*sizeof(uint32_t), cudaCpuDeviceId);

        countMin=0;
        for(int j=0; j<Ncells; j++){
            count[j] = 0xFFFFFFFFu;
        }
        // countMinStage3=0;
        for(int j=0; j<Ngrids; j++){
            countStage3[j] = 0xFFFFFFFFu;
        }

    }
    std::cout << " Stage1 KERNEL CUDA time elapsed: " << total/n_tests << " ms\n";
    
    cudaFree(count);
    cudaFree(countStage3);

	return;
}
